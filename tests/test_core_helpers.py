import numpy as np
from scipy import stats

from addivortes import _core
from addivortes.model import AddiVortesRegressor


def test_cell_indices_finds_nearest_centres_for_active_dimension():
    query = np.array([[1.0, 0.0], [9.0, 0.0], [4.9, 3.0]])
    centres = np.array([[0.0], [10.0]])
    dims = np.array([0], dtype=np.int32)
    metric_red = np.array([0], dtype=np.int32)
    member_red = np.array([2], dtype=np.int32)

    indices = _core.cell_indices(query, centres, dims, metric_red, member_red)

    np.testing.assert_array_equal(indices, [0, 1, 0])


def test_estimate_lambda_returns_positive_finite_value_near_objective_minimum():
    model = AddiVortesRegressor(nu=6.0, q=0.85)
    sigma_squared_hat = 0.2

    lambda_value = model._estimate_lambda(sigma_squared_hat)

    assert np.isfinite(lambda_value)
    assert lambda_value > 0

    def objective(candidate):
        quantile = stats.invgamma.ppf(model.q, a=model.nu / 2.0, scale=model.nu * candidate / 2.0)
        return (sigma_squared_hat - quantile) ** 2

    assert objective(lambda_value) <= objective(lambda_value * 0.5)
    assert objective(lambda_value) <= objective(lambda_value * 1.5)
