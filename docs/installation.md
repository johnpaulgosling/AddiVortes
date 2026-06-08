# Installation

AddiVortes targets Python 3.10 or newer.

## From PyPI

After a release has been published to PyPI:

```bash
python -m pip install addivortes
```

## From GitHub

Before the first PyPI release, or to install the latest development version:

```bash
python -m pip install "git+https://github.com/johnpaulgosling/AddiVortes.git"
```

To install a specific branch:

```bash
python -m pip install "git+https://github.com/johnpaulgosling/AddiVortes.git@main"
```

## From a source checkout

For local development:

```bash
git clone https://github.com/johnpaulgosling/AddiVortes.git
cd AddiVortes
python -m pip install -e ".[dev]"
```

## Build requirements

The package builds a C++20 extension with pybind11. Source installs require:

- a Python development environment with headers,
- a C++ compiler with C++20 support,
- `pip` with PEP 517 build isolation support.

Pre-built wheels avoid requiring a compiler for end users. The package checks
workflow builds and validates both wheel and source distributions.
