from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


METRIC_TO_INT = {"E": 0, "S": 1, "C": 2}


@dataclass(frozen=True)
class CategoryEncoding:
    cat_col_indices: tuple[int, ...]
    levels: dict[int, tuple[str, ...]]
    cat_scaling: float
    original_columns: tuple[str, ...]
    encoded_columns: tuple[str, ...]
    encoded_binary_cols: tuple[int, ...]


@dataclass(frozen=True)
class DesignMatrix:
    values: np.ndarray
    metric: np.ndarray
    members: np.ndarray
    columns: tuple[str, ...]
    encoding: CategoryEncoding | None


def as_dataframe(x, *, column_prefix: str = "x") -> pd.DataFrame:
    if isinstance(x, pd.DataFrame):
        return x.copy()

    arr = np.asarray(x)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError("X must be a 1- or 2-dimensional array-like object.")

    columns = [f"{column_prefix}{i}" for i in range(arr.shape[1])]
    return pd.DataFrame(arr, columns=columns)


def normalize_metric(metric, n_columns: int) -> list[str]:
    if isinstance(metric, str):
        metrics = [metric] * n_columns
    else:
        metrics = list(metric)
        if len(metrics) == 1:
            metrics = metrics * n_columns

    if len(metrics) != n_columns:
        raise ValueError("metric must be a scalar or have one value per column.")

    return [_normalize_one_metric(value) for value in metrics]


def normalize_members(members, metric: Iterable[str]) -> np.ndarray:
    metric = list(metric)
    if members is None:
        out = np.zeros(len(metric), dtype=np.int32)
        out[np.asarray(metric) == "E"] = 1
        out[np.asarray(metric) == "S"] = 2
        out[np.asarray(metric) == "C"] = 3
        return out

    arr = np.asarray(members, dtype=np.int32)
    if arr.ndim != 1 or arr.size != len(metric):
        raise ValueError("members must have one value per original covariate.")
    return arr


def prepare_design(
    x,
    *,
    metric,
    members=None,
    cat_scaling: float = 1.0,
    encoding: CategoryEncoding | None = None,
) -> DesignMatrix:
    if cat_scaling <= 0:
        raise ValueError("cat_scaling must be positive.")

    frame = as_dataframe(x)
    original_columns = tuple(str(col) for col in frame.columns)

    if encoding is not None:
        if tuple(encoding.original_columns) != original_columns:
            if len(encoding.original_columns) != frame.shape[1]:
                raise ValueError("X has a different number of columns than the training data.")
            frame.columns = list(encoding.original_columns)
            original_columns = tuple(encoding.original_columns)
        cat_scaling = encoding.cat_scaling

    metric_labels = normalize_metric(metric, frame.shape[1])
    metric_labels = [
        "C" if label == "E" and _is_categorical(frame.iloc[:, idx]) else label
        for idx, label in enumerate(metric_labels)
    ]
    member_labels = normalize_members(members, metric_labels)

    categorical_cols = _categorical_columns(frame, metric_labels, encoding)
    result_cols: list[np.ndarray] = []
    result_names: list[str] = []
    result_metric: list[int] = []
    result_members: list[int] = []
    binary_cols: list[int] = []
    levels: dict[int, tuple[str, ...]] = {}

    for col_idx, col_name in enumerate(original_columns):
        col = frame.iloc[:, col_idx]
        if col_idx in categorical_cols:
            if encoding is None:
                levels[col_idx] = _levels_for(col)
            else:
                levels[col_idx] = encoding.levels[col_idx]

            values = col.astype("string").fillna("<NA>").astype(str)
            for level in levels[col_idx][1:]:
                result_cols.append((values == level).to_numpy(dtype=float) * cat_scaling)
                result_names.append(f"{col_name}_{level}")
                result_metric.append(METRIC_TO_INT["E"])
                result_members.append(int(member_labels[col_idx]))
                binary_cols.append(len(result_cols) - 1)
        else:
            try:
                numeric = pd.to_numeric(col, errors="raise").to_numpy(dtype=float)
            except Exception as exc:  # noqa: BLE001 - provide a clearer domain error
                raise ValueError(f"Column {col_name!r} must be numeric or categorical.") from exc
            result_cols.append(numeric)
            result_names.append(col_name)
            result_metric.append(METRIC_TO_INT[metric_labels[col_idx]])
            result_members.append(int(member_labels[col_idx]))

    if not result_cols:
        raise ValueError("X must contain at least one usable covariate.")

    values = np.column_stack(result_cols).astype(float, copy=False)
    if not np.all(np.isfinite(values)):
        raise ValueError("X contains missing or non-finite values after preprocessing.")

    built_encoding = encoding
    if encoding is None and categorical_cols:
        built_encoding = CategoryEncoding(
            cat_col_indices=tuple(categorical_cols),
            levels=levels,
            cat_scaling=float(cat_scaling),
            original_columns=original_columns,
            encoded_columns=tuple(result_names),
            encoded_binary_cols=tuple(binary_cols),
        )

    return DesignMatrix(
        values=np.ascontiguousarray(values, dtype=float),
        metric=np.asarray(result_metric, dtype=np.int32),
        members=np.asarray(result_members, dtype=np.int32),
        columns=tuple(result_names),
        encoding=built_encoding,
    )


def scale_vector(y) -> tuple[np.ndarray, float, float]:
    arr = np.asarray(y, dtype=float)
    if arr.ndim != 1:
        raise ValueError("y must be a one-dimensional numeric array.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("y contains missing or non-finite values.")

    centre = float((np.max(arr) + np.min(arr)) / 2.0)
    data_range = float(np.max(arr) - np.min(arr))
    if data_range <= 0:
        raise ValueError("y has zero range and cannot be scaled.")
    return (arr - centre) / data_range, centre, data_range


def scale_matrix(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(x, dtype=float)
    centres = (np.max(values, axis=0) + np.min(values, axis=0)) / 2.0
    ranges = np.max(values, axis=0) - np.min(values, axis=0)
    if np.any(ranges <= 0):
        bad = np.where(ranges <= 0)[0].tolist()
        raise ValueError(f"X contains zero-range columns at encoded positions {bad}.")
    return (values - centres) / ranges, centres, ranges


def apply_scaling(x: np.ndarray, centres: np.ndarray, ranges: np.ndarray) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    if values.shape[1] != centres.size or values.shape[1] != ranges.size:
        raise ValueError("X has a different encoded column count than the fitted model.")
    return (values - centres) / ranges


def reduced_metric_and_members(metric: np.ndarray, members: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if metric.size != members.size:
        raise ValueError("metric and members must have the same length.")

    reduced_metric: list[int] = []
    member_counts: list[int] = []
    idx = 0
    while idx < members.size:
        current_member = members[idx]
        count = 0
        while idx + count < members.size and members[idx + count] == current_member:
            count += 1
        reduced_metric.append(int(metric[idx]))
        member_counts.append(count)
        idx += count

    return np.asarray(reduced_metric, dtype=np.int32), np.asarray(member_counts, dtype=np.int32)


def _normalize_one_metric(value) -> str:
    text = str(value).lower()
    if text in {"e", "euc", "euclidean", "0"}:
        return "E"
    if text in {"s", "sphere", "spherical", "1"}:
        return "S"
    if text in {"c", "cat", "categorical", "2"}:
        return "C"
    raise ValueError(f"Unknown metric value {value!r}.")


def _categorical_columns(
    frame: pd.DataFrame,
    metric: list[str],
    encoding: CategoryEncoding | None,
) -> list[int]:
    if encoding is not None:
        return list(encoding.cat_col_indices)

    out: list[int] = []
    for idx, label in enumerate(metric):
        col = frame.iloc[:, idx]
        if label == "C":
            if not _is_categorical(col):
                raise ValueError(
                    f"Column {frame.columns[idx]!r} is marked categorical; use object, string, "
                    "or pandas Categorical dtype."
                )
            out.append(idx)
        elif _is_categorical(col):
            out.append(idx)
    return out


def _is_categorical(series: pd.Series) -> bool:
    return (
        pd.api.types.is_object_dtype(series)
        or pd.api.types.is_string_dtype(series)
        or isinstance(series.dtype, pd.CategoricalDtype)
    )


def _levels_for(series: pd.Series) -> tuple[str, ...]:
    if isinstance(series.dtype, pd.CategoricalDtype):
        levels = [str(value) for value in series.cat.categories]
    else:
        values = series.astype("string").fillna("<NA>").astype(str)
        levels = sorted(values.unique().tolist())
    if not levels:
        raise ValueError("Categorical columns must contain at least one level.")
    return tuple(levels)
