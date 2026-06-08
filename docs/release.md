# Release checklist

This repository is configured for PyPI-oriented Python releases.

## Before release

Run locally:

```bash
python -m pip install -e ".[dev]"
python -m pytest
python -m build
python -m twine check dist/*
mkdocs build --strict
```

## GitHub checks

Pull requests run:

- tests across supported Python versions,
- wheel and source distribution builds,
- pip installation from the built wheel,
- pip installation from the built source distribution,
- strict documentation builds.

## Publishing to PyPI

The `Publish Python package` workflow runs on GitHub releases and uses PyPI
trusted publishing.

Repository maintainers need to configure a PyPI trusted publisher for this
repository and the `pypi` GitHub environment before the first release.

After the release workflow succeeds, users can install with:

```bash
python -m pip install addivortes
```
