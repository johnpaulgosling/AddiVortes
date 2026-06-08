import numpy as np
import pandas as pd
import pytest

from addivortes import AddiVortes, AddiVortesRegressor

from conftest import fast_model, regression_data


def test_import_aliases_are_available():
    assert AddiVortes is AddiVortesRegressor


def test_numeric_fit_predict_response_and_quantiles():
    x, y = regression_data()

    model = fast_model().fit(x, y)

    pred = model.predict(x[:5])
    assert pred.shape == (5,)
    assert np.all(np.isfinite(pred))

    quantiles = model.predict(x[:5], kind="quantile", quantiles=(0.1, 0.9))
    assert quantiles.shape == (5, 2)
    assert np.all(np.isfinite(quantiles))
    assert model.summary()["posterior_samples"] == 5


def test_summary_repr_score_and_fit_predict():
    x, y = regression_data(seed=321)
    model = fast_model().fit(x, y)

    summary = model.summary()

    assert "AddiVortesRegressor" in repr(model)
    assert summary["n_features"] == 3
    assert summary["n_tessellations"] == 4
    assert summary["posterior_samples"] == 5
    assert np.isfinite(summary["in_sample_rmse"])
    assert np.isfinite(summary["lambda"])
    assert np.isfinite(model.score(x, y))

    fit_predict = fast_model(random_state=123).fit_predict(x, y)
    assert fit_predict.shape == (x.shape[0],)
    assert np.all(np.isfinite(fit_predict))


def test_categorical_dataframe_unknown_level_prediction():
    rng = np.random.default_rng(456)
    frame = pd.DataFrame(
        {
            "x1": rng.normal(size=24),
            "group": pd.Categorical(rng.choice(["a", "b", "c"], size=24), categories=["a", "b", "c"]),
        }
    )
    y = frame["x1"].to_numpy() + (frame["group"].astype(str) == "b").to_numpy(dtype=float)

    model = AddiVortesRegressor(
        n_tessellations=3,
        total_mcmc_iter=14,
        burn_in=4,
        thinning=2,
        random_state=7,
    ).fit(frame, y)

    new_frame = pd.DataFrame(
        {
            "x1": [0.0, 1.0],
            "group": pd.Categorical(["b", "d"], categories=["a", "b", "c", "d"]),
        }
    )
    pred = model.predict(new_frame)
    assert pred.shape == (2,)
    assert np.all(np.isfinite(pred))


def test_sklearn_style_params_round_trip():
    model = AddiVortesRegressor(n_tessellations=2)
    assert model.get_params()["n_tessellations"] == 2
    model.set_params(n_tessellations=3)
    assert model.n_tessellations == 3

    with pytest.raises(ValueError, match="Invalid parameter"):
        model.set_params(not_a_parameter=True)


def test_predict_validation_errors():
    x, y = regression_data()
    model = fast_model().fit(x, y)

    with pytest.raises(RuntimeError, match="not fitted"):
        fast_model().predict(x)

    with pytest.raises(ValueError, match="different number of columns"):
        model.predict(x[:, :2])

    with pytest.raises(ValueError, match="kind"):
        model.predict(x, kind="bad")

    with pytest.raises(ValueError, match="interval"):
        model.predict(x, interval="bad")

    with pytest.raises(ValueError, match="quantiles"):
        model.predict(x, kind="quantile", quantiles=(-0.1, 1.1))


def test_fit_validation_errors_and_high_dimension_warning():
    x, y = regression_data(n_obs=8, n_features=3)

    with pytest.raises(ValueError, match="same number of observations"):
        fast_model().fit(x, y[:-1])

    with pytest.raises(ValueError, match="burn_in"):
        AddiVortesRegressor(total_mcmc_iter=5, burn_in=5).fit(x, y)

    x_high_dim, y_high_dim = regression_data(seed=99, n_obs=5, n_features=7)
    with pytest.warns(RuntimeWarning, match="exceeds number of observations"):
        fast_model(n_tessellations=2, total_mcmc_iter=8, burn_in=4, omega=1).fit(x_high_dim, y_high_dim)
