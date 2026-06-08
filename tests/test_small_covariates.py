import numpy as np
import pytest

from conftest import fast_model, regression_data


@pytest.mark.parametrize("n_features", [1, 2, 3])
def test_small_covariate_counts_fit_and_predict(n_features):
    x, y = regression_data(seed=100 + n_features, n_obs=18, n_features=n_features)

    model = fast_model(n_tessellations=3, total_mcmc_iter=14, burn_in=4, omega=1).fit(x, y)
    pred = model.predict(x[:5])

    assert np.isfinite(model.in_sample_rmse_)
    assert model.in_sample_rmse_ >= 0
    assert model.posterior_.sigma.shape == (5,)
    assert np.all(np.isfinite(model.posterior_.sigma))
    assert np.all(model.posterior_.sigma > 0)
    assert pred.shape == (5,)
    assert np.all(np.isfinite(pred))
