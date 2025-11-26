# CombAero Python bindings

This package provides Python bindings for the core CombAero C++ library:

- Thermodynamic and transport properties for gas mixtures
- Combustion helpers (oxygen requirement, complete combustion)
- Equivalence ratio and Bilger mixture-fraction utilities
- Humid air properties and standard dry-air composition
- Species common-name mapping

The bindings are implemented with [pybind11](https://pybind11.readthedocs.io/) and built via [scikit-build-core](https://scikit-build-core.readthedocs.io/).

## Installation (from source)

The recommended workflow is from the repository root; see the top-level README for details. In short:

```bash
python -m build --wheel
python -m pip install dist/combaero-0.0.1-*.whl
```

## Usage

Example using NumPy arrays and the high-level `combaero` package:

```python
import numpy as np
import combaero as ca

T = 300.0
P = 101325.0
X = np.zeros(len(ca.species_common_names()), dtype=float)
X[ca.species_common_names().keys().__iter__().__next__()]  # set your composition here

print("cp =", ca.cp(T, X))
print("rho =", ca.density(T, P, X))
print("species names =", ca.species_common_names())
```
