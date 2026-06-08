import numpy as np
import pandas as pd
import pytest

from addivortes.preprocessing import (
    apply_scaling,
    normalize_metric,
    prepare_design,
    reduced_metric_and_members,
    scale_matrix,
    scale_vector,
)


def test_scale_vector_midpoint_range_scaling():
    scaled, centre, data_range = scale_vector(np.array([0.0, 5.0, 10.0]))

    assert centre == 5.0
    assert data_range == 10.0
    np.testing.assert_allclose(scaled, [-0.5, 0.0, 0.5])


def test_scale_vector_rejects_zero_range_and_non_finite_values():
    with pytest.raises(ValueError, match="zero range"):
        scale_vector(np.array([2.0, 2.0]))

    with pytest.raises(ValueError, match="non-finite"):
        scale_vector(np.array([1.0, np.nan]))


def test_scale_matrix_and_apply_scaling():
    x = np.array([[0.0, 10.0], [5.0, 20.0], [10.0, 30.0]])

    scaled, centres, ranges = scale_matrix(x)

    np.testing.assert_allclose(centres, [5.0, 20.0])
    np.testing.assert_allclose(ranges, [10.0, 20.0])
    np.testing.assert_allclose(scaled, [[-0.5, -0.5], [0.0, 0.0], [0.5, 0.5]])
    np.testing.assert_allclose(apply_scaling(x, centres, ranges), scaled)


def test_scale_matrix_and_apply_scaling_validation():
    with pytest.raises(ValueError, match="zero-range"):
        scale_matrix(np.array([[1.0, 2.0], [1.0, 3.0]]))

    with pytest.raises(ValueError, match="different encoded column count"):
        apply_scaling(np.ones((2, 3)), np.ones(2), np.ones(2))


def test_prepare_design_encodes_character_categorical_columns():
    frame = pd.DataFrame({"x1": [0.0, 1.0, 2.0], "group": ["b", "a", "c"]})

    design = prepare_design(frame, metric="euclidean", cat_scaling=2.0)

    assert design.columns == ("x1", "group_b", "group_c")
    assert design.encoding is not None
    assert design.encoding.cat_col_indices == (1,)
    assert design.encoding.levels[1] == ("a", "b", "c")
    assert design.encoding.encoded_binary_cols == (1, 2)
    np.testing.assert_allclose(design.values, [[0.0, 2.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 2.0]])


def test_prepare_design_reuses_encoding_and_maps_unseen_levels_to_reference():
    train = pd.DataFrame({"x1": [0.0, 1.0, 2.0], "group": ["a", "b", "c"]})
    design = prepare_design(train, metric="euclidean", cat_scaling=1.5)
    new = pd.DataFrame({"x1": [3.0, 4.0], "group": ["b", "unknown"]})

    encoded = prepare_design(new, metric="euclidean", encoding=design.encoding)

    assert encoded.columns == ("x1", "group_b", "group_c")
    np.testing.assert_allclose(encoded.values, [[3.0, 1.5, 0.0], [4.0, 0.0, 0.0]])


def test_prepare_design_supports_multiple_and_categorical_only_columns():
    frame = pd.DataFrame(
        {
            "group": pd.Categorical(["low", "mid", "high"], categories=["low", "mid", "high"]),
            "flag": ["yes", "no", "yes"],
        }
    )

    design = prepare_design(frame, metric="categorical")

    assert design.columns == ("group_mid", "group_high", "flag_yes")
    assert design.encoding is not None
    assert design.values.shape == (3, 3)
    assert set(np.unique(design.values)) == {0.0, 1.0}


def test_prepare_design_validation():
    with pytest.raises(ValueError, match="cat_scaling"):
        prepare_design([[1.0], [2.0]], metric="euclidean", cat_scaling=0.0)

    with pytest.raises(ValueError, match="marked categorical"):
        prepare_design(pd.DataFrame({"x": [1.0, 2.0]}), metric="categorical")

    with pytest.raises(ValueError, match="Unknown metric"):
        normalize_metric("bad-metric", 1)


def test_reduced_metric_and_members_groups_contiguous_memberships():
    metric = np.array([0, 0, 1, 1, 0], dtype=np.int32)
    members = np.array([1, 1, 2, 2, 3], dtype=np.int32)

    reduced_metric, member_counts = reduced_metric_and_members(metric, members)

    np.testing.assert_array_equal(reduced_metric, [0, 1, 0])
    np.testing.assert_array_equal(member_counts, [2, 2, 1])
