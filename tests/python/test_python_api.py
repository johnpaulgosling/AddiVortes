import numpy as np
import pandas as pd

from addivortes import AddiVortesRegressor


def test_numeric_fit_predict_response_and_quantiles():
    rng = np.random.default_rng(123)
    x = rng.normal(size=(18, 3))
    y = x[:, 0] - 0.5 * x[:, 1] + rng.normal(scale=0.1, size=18)

    model = AddiVortesRegressor(
        n_tessellations=4,
        total_mcmc_iter=16,
        burn_in=6,
        thinning=2,
        random_state=42,
    ).fit(x, y)

    pred = model.predict(x[:5])
    assert pred.shape == (5,)
    assert np.all(np.isfinite(pred))

    quantiles = model.predict(x[:5], kind="quantile", quantiles=(0.1, 0.9))
    assert quantiles.shape == (5, 2)
    assert np.all(np.isfinite(quantiles))
    assert model.summary()["posterior_samples"] == 5


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
