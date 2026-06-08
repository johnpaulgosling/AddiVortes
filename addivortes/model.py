from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import warnings

import numpy as np
from scipy import optimize, stats

from . import _core
from .preprocessing import (
    CategoryEncoding,
    apply_scaling,
    normalize_members,
    normalize_metric,
    prepare_design,
    reduced_metric_and_members,
    scale_matrix,
    scale_vector,
)


@dataclass(frozen=True)
class PosteriorSamples:
    tessellations: list[list[np.ndarray]]
    dimensions: list[list[np.ndarray]]
    predictions: list[list[np.ndarray]]
    sigma: np.ndarray
    prediction_matrix: np.ndarray


class AddiVortesRegressor:
    """Bayesian Additive Voronoi Tessellation regressor.

    The estimator follows the usual Python machine-learning pattern:
    instantiate with model hyperparameters, call :meth:`fit`, then call
    :meth:`predict`.
    """

    def __init__(
        self,
        *,
        n_tessellations: int = 200,
        total_mcmc_iter: int = 1200,
        burn_in: int = 200,
        thinning: int = 1,
        nu: float = 6.0,
        q: float = 0.85,
        k: int = 3,
        omega: float | None = None,
        lambda_rate: float = 25.0,
        initial_sigma: str = "linear",
        metric: str | list[str] = "euclidean",
        members: list[int] | np.ndarray | None = None,
        cat_scaling: float = 1.0,
        random_state: int | None = None,
        verbose: bool = False,
    ) -> None:
        self.n_tessellations = n_tessellations
        self.total_mcmc_iter = total_mcmc_iter
        self.burn_in = burn_in
        self.thinning = thinning
        self.nu = nu
        self.q = q
        self.k = k
        self.omega = omega
        self.lambda_rate = lambda_rate
        self.initial_sigma = initial_sigma
        self.metric = metric
        self.members = members
        self.cat_scaling = cat_scaling
        self.random_state = random_state
        self.verbose = verbose

    def fit(self, X, y) -> "AddiVortesRegressor":
        """Fit the AddiVortes regression model."""

        self._validate_hyperparameters()

        raw_columns = np.asarray(X).shape[1] if np.asarray(X).ndim == 2 else 1
        metric_labels = normalize_metric(self.metric, raw_columns)
        original_members = normalize_members(self.members, metric_labels)
        design = prepare_design(
            X,
            metric=metric_labels,
            members=original_members,
            cat_scaling=float(self.cat_scaling),
        )
        y_scaled, y_centre, y_range = scale_vector(y)

        if design.values.shape[0] != y_scaled.size:
            raise ValueError("X and y must contain the same number of observations.")

        x_scaled, x_centres, x_ranges = scale_matrix(design.values)
        x_scaled = self._restore_unscaled_columns(x_scaled, design.values, design.metric, design.encoding)

        n_obs, n_features = x_scaled.shape
        if n_features > n_obs:
            warnings.warn(
                "Number of covariates exceeds number of observations; model results may not be stable.",
                RuntimeWarning,
                stacklevel=2,
            )
        omega = min(3.0, float(raw_columns)) if self.omega is None else float(self.omega)
        sigma_squared_mu = (0.5 / (float(self.k) * np.sqrt(float(self.n_tessellations)))) ** 2
        sigma_squared_hat = self._initial_sigma_squared(x_scaled, y_scaled)
        lambda_value = self._estimate_lambda(sigma_squared_hat)

        proposal_sd = self._proposal_sd(x_ranges, design.metric)
        proposal_mu = np.zeros(n_features, dtype=float)
        proposal_mu[design.metric != 0] = x_centres[design.metric != 0]

        rng = np.random.default_rng(self.random_state)
        init_tess, init_dim, init_pred = self._initial_state(
            rng,
            n_features=n_features,
            y_scaled=y_scaled,
            proposal_sd=proposal_sd,
            metric=design.metric,
            members=design.members,
            encoding=design.encoding,
        )

        binary_cols = (
            np.asarray(design.encoding.encoded_binary_cols, dtype=np.int32)
            if design.encoding is not None
            else np.empty(0, dtype=np.int32)
        )
        seed = int(rng.integers(0, np.iinfo(np.uint64).max, dtype=np.uint64))

        result = _core.run_mcmc(
            x_scaled,
            y_scaled,
            design.metric,
            design.members,
            int(self.n_tessellations),
            int(self.total_mcmc_iter),
            int(self.burn_in),
            int(self.thinning),
            float(self.nu),
            float(lambda_value),
            float(sigma_squared_mu),
            float(omega),
            float(self.lambda_rate),
            proposal_sd,
            proposal_mu,
            init_tess,
            init_dim,
            init_pred,
            binary_cols,
            float(self.cat_scaling),
            seed,
            bool(self.verbose),
        )

        posterior = PosteriorSamples(
            tessellations=result["posterior_tess"],
            dimensions=result["posterior_dim"],
            predictions=result["posterior_pred"],
            sigma=np.asarray(result["posterior_sigma"], dtype=float),
            prediction_matrix=np.asarray(result["prediction_matrix"], dtype=float),
        )

        if posterior.prediction_matrix.shape[1] > 0:
            mean_yhat = posterior.prediction_matrix.mean(axis=1) * y_range + y_centre
            y_arr = np.asarray(y, dtype=float)
            in_sample_rmse = float(np.sqrt(np.mean((y_arr - mean_yhat) ** 2)))
        else:
            in_sample_rmse = float("nan")

        self.posterior_ = posterior
        self.x_centres_ = x_centres
        self.x_ranges_ = x_ranges
        self.y_centre_ = y_centre
        self.y_range_ = y_range
        self.in_sample_rmse_ = in_sample_rmse
        self.metric_ = design.metric
        self.members_ = design.members
        self.metric_labels_ = metric_labels
        self.original_members_ = original_members
        self.cat_encoding_: CategoryEncoding | None = design.encoding
        self.feature_names_in_ = tuple(design.columns)
        self.n_original_features_in_ = raw_columns
        self.n_features_in_ = n_features
        self.metric_red_, self.member_red_ = reduced_metric_and_members(design.metric, design.members)
        self.lambda_ = float(lambda_value)
        self.sigma_squared_mu_ = float(sigma_squared_mu)
        return self

    def predict(
        self,
        X,
        *,
        kind: str = "response",
        quantiles: tuple[float, ...] | list[float] = (0.025, 0.975),
        interval: str = "credible",
        random_state: int | None = None,
    ) -> np.ndarray:
        """Predict responses or posterior quantiles for new observations."""

        self._check_is_fitted()
        kind = _match_choice(kind, {"response", "quantile"}, "kind")
        interval = _match_choice(interval, {"credible", "prediction"}, "interval")
        quantiles_arr = np.asarray(quantiles, dtype=float)
        if np.any((quantiles_arr < 0) | (quantiles_arr > 1)):
            raise ValueError("quantiles must be probabilities in [0, 1].")

        raw_columns = np.asarray(X).shape[1] if np.asarray(X).ndim == 2 else 1
        if raw_columns != self.n_original_features_in_:
            raise ValueError("X has a different number of columns than the fitted model.")

        design = prepare_design(
            X,
            metric=self.metric_labels_,
            members=self.original_members_,
            cat_scaling=float(self.cat_scaling),
            encoding=self.cat_encoding_,
        )
        if design.values.shape[1] != self.n_features_in_:
            raise ValueError("X has a different encoded feature count than the fitted model.")

        x_scaled = apply_scaling(design.values, self.x_centres_, self.x_ranges_)
        x_scaled = self._restore_unscaled_columns(x_scaled, design.values, self.metric_, self.cat_encoding_)

        n_obs = x_scaled.shape[0]
        n_samples = len(self.posterior_.tessellations)
        if n_samples == 0:
            raise RuntimeError("The fitted model contains no posterior samples.")

        rng = np.random.default_rng(self.random_state if random_state is None else random_state)
        prediction_matrix = np.empty((n_obs, n_samples), dtype=float)

        for sample_idx in range(n_samples):
            sample_pred = np.zeros(n_obs, dtype=float)
            tess_sample = self.posterior_.tessellations[sample_idx]
            dim_sample = self.posterior_.dimensions[sample_idx]
            pred_sample = self.posterior_.predictions[sample_idx]

            for tess, dims, cell_pred in zip(tess_sample, dim_sample, pred_sample, strict=True):
                indices = _core.cell_indices(
                    x_scaled,
                    np.asarray(tess, dtype=float),
                    np.asarray(dims, dtype=np.int32),
                    self.metric_red_,
                    self.member_red_,
                )
                sample_pred += np.asarray(cell_pred, dtype=float)[indices]

            if kind == "quantile" and interval == "prediction":
                sample_pred += rng.normal(
                    loc=0.0,
                    scale=float(np.sqrt(self.posterior_.sigma[sample_idx])),
                    size=n_obs,
                )

            prediction_matrix[:, sample_idx] = sample_pred

        if kind == "response":
            return prediction_matrix.mean(axis=1) * self.y_range_ + self.y_centre_

        return np.quantile(prediction_matrix, quantiles_arr, axis=1).T * self.y_range_ + self.y_centre_

    def fit_predict(self, X, y, **predict_kwargs) -> np.ndarray:
        """Fit the model and return predictions for the training data."""

        return self.fit(X, y).predict(X, **predict_kwargs)

    def score(self, X, y) -> float:
        """Return the coefficient of determination on test data."""

        y_arr = np.asarray(y, dtype=float)
        pred = self.predict(X)
        ss_res = float(np.sum((y_arr - pred) ** 2))
        ss_tot = float(np.sum((y_arr - np.mean(y_arr)) ** 2))
        return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    def summary(self) -> dict[str, Any]:
        """Return a dictionary of fitted model summary information."""

        self._check_is_fitted()
        tess_sizes = [
            [np.asarray(tess).shape[0] for tess in sample]
            for sample in self.posterior_.tessellations[: min(5, len(self.posterior_.tessellations))]
        ]
        flat_sizes = [size for sample in tess_sizes for size in sample]
        return {
            "n_features": self.n_features_in_,
            "n_tessellations": self.n_tessellations,
            "posterior_samples": len(self.posterior_.tessellations),
            "in_sample_rmse": self.in_sample_rmse_,
            "lambda": self.lambda_,
            "average_cells_per_tessellation": float(np.mean(flat_sizes)) if flat_sizes else float("nan"),
        }

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        del deep
        return {
            "n_tessellations": self.n_tessellations,
            "total_mcmc_iter": self.total_mcmc_iter,
            "burn_in": self.burn_in,
            "thinning": self.thinning,
            "nu": self.nu,
            "q": self.q,
            "k": self.k,
            "omega": self.omega,
            "lambda_rate": self.lambda_rate,
            "initial_sigma": self.initial_sigma,
            "metric": self.metric,
            "members": self.members,
            "cat_scaling": self.cat_scaling,
            "random_state": self.random_state,
            "verbose": self.verbose,
        }

    def set_params(self, **params) -> "AddiVortesRegressor":
        valid = self.get_params()
        for key, value in params.items():
            if key not in valid:
                raise ValueError(f"Invalid parameter {key!r}.")
            setattr(self, key, value)
        return self

    def __repr__(self) -> str:
        params = self.get_params()
        compact = ", ".join(f"{key}={value!r}" for key, value in params.items())
        return f"AddiVortesRegressor({compact})"

    def _validate_hyperparameters(self) -> None:
        if int(self.n_tessellations) <= 0:
            raise ValueError("n_tessellations must be positive.")
        if int(self.total_mcmc_iter) <= 0:
            raise ValueError("total_mcmc_iter must be positive.")
        if int(self.burn_in) < 0:
            raise ValueError("burn_in must be non-negative.")
        if int(self.thinning) <= 0:
            raise ValueError("thinning must be positive.")
        if int(self.burn_in) >= int(self.total_mcmc_iter):
            raise ValueError("burn_in must be smaller than total_mcmc_iter.")
        if float(self.nu) <= 0:
            raise ValueError("nu must be positive.")
        if not (0 < float(self.q) < 1):
            raise ValueError("q must be in (0, 1).")
        if int(self.k) <= 0:
            raise ValueError("k must be positive.")
        if float(self.lambda_rate) <= 0:
            raise ValueError("lambda_rate must be positive.")
        if float(self.cat_scaling) <= 0:
            raise ValueError("cat_scaling must be positive.")

    def _initial_sigma_squared(self, x_scaled: np.ndarray, y_scaled: np.ndarray) -> float:
        if str(self.initial_sigma).lower() == "naive":
            return float(np.var(y_scaled, ddof=1))

        design = np.column_stack([np.ones(x_scaled.shape[0]), x_scaled])
        coef, *_ = np.linalg.lstsq(design, y_scaled, rcond=None)
        resid = y_scaled - design @ coef
        dof = max(1, x_scaled.shape[0] - x_scaled.shape[1] - 1)
        return float(np.sum(resid**2) / dof)

    def _estimate_lambda(self, sigma_squared_hat: float) -> float:
        shape = float(self.nu) / 2.0

        def objective(lambda_value: float) -> float:
            rate = float(self.nu) * lambda_value / 2.0
            quantile = stats.invgamma.ppf(float(self.q), a=shape, scale=rate)
            return float((sigma_squared_hat - quantile) ** 2)

        result = optimize.minimize_scalar(objective, bounds=(0.001, 100.0), method="bounded")
        return float(result.x)

    @staticmethod
    def _proposal_sd(ranges: np.ndarray, metric: np.ndarray) -> np.ndarray:
        sd = (ranges / 2.0) / stats.norm.ppf(0.75)
        sd = np.asarray(sd, dtype=float)
        sd[metric == 0] = 0.8
        return sd

    def _initial_state(
        self,
        rng: np.random.Generator,
        *,
        n_features: int,
        y_scaled: np.ndarray,
        proposal_sd: np.ndarray,
        metric: np.ndarray,
        members: np.ndarray,
        encoding: CategoryEncoding | None,
    ) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
        init_tess: list[np.ndarray] = []
        init_dim: list[np.ndarray] = []
        init_pred: list[np.ndarray] = []
        binary_cols = set(encoding.encoded_binary_cols if encoding is not None else ())

        for _ in range(int(self.n_tessellations)):
            dim = np.asarray([int(rng.integers(0, n_features))], dtype=np.int32)
            col = dim[0]
            value = rng.normal(loc=0.0, scale=proposal_sd[col])
            if metric[col] == 1 and _is_last_member_column(col, members):
                value = _period_shift(value, np.pi)
            if col in binary_cols:
                value = rng.uniform(0.0, float(self.cat_scaling))

            init_dim.append(dim)
            init_tess.append(np.asarray([[value]], dtype=float))
            init_pred.append(np.asarray([float(np.mean(y_scaled) / float(self.n_tessellations))], dtype=float))

        return init_tess, init_dim, init_pred

    @staticmethod
    def _restore_unscaled_columns(
        scaled: np.ndarray,
        original: np.ndarray,
        metric: np.ndarray,
        encoding: CategoryEncoding | None,
    ) -> np.ndarray:
        out = np.asarray(scaled, dtype=float).copy()
        out[:, metric != 0] = original[:, metric != 0]
        if encoding is not None and encoding.encoded_binary_cols:
            binary_cols = np.asarray(encoding.encoded_binary_cols, dtype=int)
            out[:, binary_cols] = original[:, binary_cols]
        return np.ascontiguousarray(out, dtype=float)

    def _check_is_fitted(self) -> None:
        if not hasattr(self, "posterior_"):
            raise RuntimeError("This AddiVortesRegressor instance is not fitted yet.")


def _period_shift(value: float, limit: float) -> float:
    while value >= limit:
        value -= 2 * limit
    while value < -limit:
        value += 2 * limit
    return value


def _is_last_member_column(index: int, members: np.ndarray) -> bool:
    return index == members.size - 1 or members[index + 1] != members[index]


def _match_choice(value: str, choices: set[str], name: str) -> str:
    normalized = str(value).lower()
    if normalized not in choices:
        raise ValueError(f"{name} must be one of {sorted(choices)}.")
    return normalized


AddiVortes = AddiVortesRegressor
