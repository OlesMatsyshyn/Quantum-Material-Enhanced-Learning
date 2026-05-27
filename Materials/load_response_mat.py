from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat


REQUIRED_FIELDS = (
    "material_name",
    "model_name",
    "version",
    "mu_grid",
    "S_photo",
    "Omega_z_photo",
    "D_photo",
    "J_photo",
    "units",
    "parameters",
    "grid",
    "index_convention",
)

REQUIRED_UNITS = {
    "energy": "meV",
    "length": "nm",
    "momentum": "nm^{-1}",
    "S": "meV",
    "Omega_z": "dimensionless",
    "D": "nm",
    "J": "meV nm",
}


def _as_1d(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=float).reshape(-1)


def _interp_array(mu_grid: np.ndarray, values: np.ndarray, mu: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    flat = values.reshape((values.shape[0], -1))
    out = np.array([np.interp(mu, mu_grid, flat[:, i]) for i in range(flat.shape[1])])
    return out.reshape(values.shape[1:])


def validate_response_data(data: dict[str, Any]) -> None:
    missing = [field for field in REQUIRED_FIELDS if field not in data]
    if missing:
        raise ValueError(f"Missing required response fields: {missing}")

    if data["version"] != "v1_response_objects":
        raise ValueError("version must be 'v1_response_objects'")

    units = data["units"]
    for key, expected in REQUIRED_UNITS.items():
        actual = units.get(key) if isinstance(units, dict) else getattr(units, key, None)
        if actual != expected:
            raise ValueError(f"units.{key} must be {expected!r}, got {actual!r}")

    mu_grid = _as_1d(data["mu_grid"])
    nmu = mu_grid.size

    if np.asarray(data["S_photo"]).shape != (nmu, 2, 2):
        raise ValueError("S_photo must have shape [Nmu, 2, 2]")
    if np.asarray(data["Omega_z_photo"]).reshape(-1).shape != (nmu,):
        raise ValueError("Omega_z_photo must have shape [Nmu] or [Nmu, 1]")
    if np.asarray(data["D_photo"]).shape != (nmu, 2):
        raise ValueError("D_photo must have shape [Nmu, 2]")
    if np.asarray(data["J_photo"]).shape != (nmu, 2, 2, 2):
        raise ValueError("J_photo must have shape [Nmu, 2, 2, 2]")

    definitions = data.get("definitions", {})
    d_definition = definitions.get("D", "") if isinstance(definitions, dict) else ""
    index_convention = str(data.get("index_convention", ""))
    if "- sum_n int_k Omega_n" not in d_definition and "- int Omega" not in index_convention:
        raise ValueError("D_photo sign convention must be documented in definitions.D or index_convention")


def load_response_mat(path: str | Path, mu: float | None = None) -> dict[str, Any]:
    data = loadmat(path, simplify_cells=True)
    data = {key: value for key, value in data.items() if not key.startswith("__")}
    validate_response_data(data)

    if mu is None:
        return data

    mu_grid = _as_1d(data["mu_grid"])
    return {
        "S": _interp_array(mu_grid, data["S_photo"], mu),
        "D": _interp_array(mu_grid, data["D_photo"], mu),
        "J": _interp_array(mu_grid, data["J_photo"], mu),
        "Omega_z": float(np.interp(mu, mu_grid, _as_1d(data["Omega_z_photo"]))),
        "parameters": data["parameters"],
        "units": data["units"],
        "material_name": data["material_name"],
        "model_name": data["model_name"],
        "mu": mu,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Load and validate a material response-object .mat file.")
    parser.add_argument("path", type=Path)
    parser.add_argument("--mu", type=float, default=None)
    args = parser.parse_args()

    loaded = load_response_mat(args.path, args.mu)
    print(f"Loaded {loaded.get('material_name')} / {loaded.get('model_name')}")
    print(f"Units: {loaded.get('units')}")
    if args.mu is not None:
        print(f"Interpolated at mu = {args.mu} meV")
        print(f"S shape: {loaded['S'].shape}, D shape: {loaded['D'].shape}, J shape: {loaded['J'].shape}")
