import numpy as np
import pandas as pd
import pytest

from conftest import fast_model


def test_model_stores_categorical_encoding_and_predicts_unseen_levels():
    rng = np.random.default_rng(456)
    frame = pd.DataFrame(
        {
            "x1": rng.normal(size=24),
            "group": pd.Categorical(rng.choice(["a", "b", "c"], size=24), categories=["a", "b", "c"]),
            "flag": rng.choice(["yes", "no"], size=24),
        }
    )
    y = frame["x1"].to_numpy() + (frame["group"].astype(str) == "b").to_numpy(dtype=float)

    model = fast_model(n_tessellations=3, total_mcmc_iter=14, burn_in=4, cat_scaling=1.75).fit(frame, y)

    assert model.cat_encoding_ is not None
    assert model.cat_encoding_.cat_scaling == 1.75
    assert model.cat_encoding_.original_columns == ("x1", "group", "flag")
    assert model.cat_encoding_.encoded_binary_cols
    assert any(name.startswith("group_") for name in model.feature_names_in_)

    new_frame = pd.DataFrame(
        {
            "x1": [0.0, 1.0],
            "group": pd.Categorical(["b", "unknown"], categories=["a", "b", "c", "unknown"]),
            "flag": ["yes", "maybe"],
        }
    )
    pred = model.predict(new_frame)

    assert pred.shape == (2,)
    assert np.all(np.isfinite(pred))


def test_binary_categorical_tessellation_centres_are_clamped_to_cat_scaling():
    frame = pd.DataFrame(
        {
            "group": ["a", "b", "c"] * 10,
        }
    )
    y = (frame["group"] == "c").to_numpy(dtype=float)
    cat_scaling = 0.6

    model = fast_model(
        n_tessellations=5,
        total_mcmc_iter=18,
        burn_in=6,
        cat_scaling=cat_scaling,
        random_state=11,
    ).fit(frame, y)

    binary_cols = set(model.cat_encoding_.encoded_binary_cols)
    observed_binary_dim = False
    for sample_tess, sample_dims in zip(model.posterior_.tessellations, model.posterior_.dimensions, strict=True):
        for tess, dims in zip(sample_tess, sample_dims, strict=True):
            for local_idx, global_dim in enumerate(dims):
                if int(global_dim) not in binary_cols:
                    continue
                observed_binary_dim = True
                values = np.asarray(tess)[:, local_idx]
                assert np.all(values >= 0.0)
                assert np.all(values <= cat_scaling)

    assert observed_binary_dim


def test_invalid_cat_scaling_rejected():
    frame = pd.DataFrame({"x1": [0.0, 1.0, 2.0], "group": ["a", "b", "a"]})
    y = np.array([0.0, 1.0, 0.5])

    with pytest.raises(ValueError, match="cat_scaling"):
        fast_model(cat_scaling=-1).fit(frame, y)
