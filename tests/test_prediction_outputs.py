import numpy as np

from conftest import fast_model, regression_data


def test_quantile_predictions_preserve_requested_order():
    x, y = regression_data(seed=777, n_obs=22, n_features=3)
    model = fast_model(n_tessellations=5, total_mcmc_iter=20, burn_in=6, random_state=777).fit(x, y)

    increasing = model.predict(x[:4], kind="quantile", quantiles=(0.25, 0.5, 0.75))
    decreasing = model.predict(x[:4], kind="quantile", quantiles=(0.95, 0.75))

    assert increasing.shape == (4, 3)
    assert decreasing.shape == (4, 2)
    assert np.all(increasing[:, 0] <= increasing[:, 1])
    assert np.all(increasing[:, 1] <= increasing[:, 2])
    assert np.all(decreasing[:, 0] >= decreasing[:, 1])


def test_credible_interval_is_default_quantile_interval():
    x, y = regression_data(seed=888, n_obs=22, n_features=3)
    model = fast_model(n_tessellations=5, total_mcmc_iter=20, burn_in=6, random_state=888).fit(x, y)

    default_interval = model.predict(x[:6], kind="quantile", quantiles=(0.1, 0.9))
    explicit_credible = model.predict(x[:6], kind="quantile", interval="credible", quantiles=(0.1, 0.9))

    np.testing.assert_allclose(default_interval, explicit_credible)


def test_prediction_intervals_include_sigma_noise_and_response_ignores_interval():
    x, y = regression_data(seed=999, n_obs=24, n_features=3)
    model = fast_model(n_tessellations=5, total_mcmc_iter=22, burn_in=6, random_state=999).fit(x, y)

    assert model.posterior_.sigma.shape == (8,)
    assert np.all(np.isfinite(model.posterior_.sigma))
    assert np.all(model.posterior_.sigma > 0)

    response = model.predict(x[:6])
    response_with_interval = model.predict(x[:6], interval="prediction")
    np.testing.assert_allclose(response, response_with_interval)

    credible = model.predict(x[:6], kind="quantile", interval="credible", quantiles=(0.1, 0.9))

    # Inflate sigma to make the effect deterministic and isolate prediction
    # interval behavior from sampler variability.
    model.posterior_.sigma[:] = 4.0
    prediction = model.predict(
        x[:6],
        kind="quantile",
        interval="prediction",
        quantiles=(0.1, 0.9),
        random_state=123,
    )

    credible_width = np.mean(credible[:, 1] - credible[:, 0])
    prediction_width = np.mean(prediction[:, 1] - prediction[:, 0])
    assert prediction_width > credible_width
