"""
Thin-film harmonic transmission from a standardized material response file.

Example command-line use:

    python transmission_from_response.py MBT/MBT_response_objects.mat \
        --omega-meV 10 \
        --mu-meV 0 \
        --gamma-meV 1 \
        --E0-Vm 100000 \
        --pol 1 0

Inputs:
    response_mat
        Path to a standardized response-object .mat file.

    --omega-meV
        Photon energy hbar*omega in meV.

    --mu-meV
        Chemical potential in meV. The code interpolates the stored response
        tensors to this value.

    --gamma-meV
        Relaxation/broadening scale Gamma in meV.

    --E0-Vm
        Incident electric-field amplitude in V/m.

    --pol px py
        Incident polarization vector. The code normalizes it internally.
        Examples:
            --pol 1 0   means x-polarized light
            --pol 0 1   means y-polarized light
            --pol 1 1   means 45-degree linear polarization

Typical material files:

    MBT:
        MBT/MBT_response_objects.mat

    WTe2:
        WTe2/WTe2_ShiSong_6band_deltaZ_response_objects.mat

    strained bilayer graphene:
        Bilayer graphene/simple_strained_BLG_photocurrent_objects_AB_xi1.mat

    twisted bilayer graphene:
        Twisted bilayer graphene (TBG)/tbg_photocurrent_objects.mat

Expected response-object fields:
    S        or S_photo          linear Drude/stiffness tensor, meV
    Omega_z  or Omega_z_photo    integrated Berry curvature, dimensionless
    D        or D_photo          Berry-curvature dipole vector, nm
    J        or J_photo          jerk tensor, meV nm
    mu_grid                      chemical-potential grid, meV

Physics included:
    linear response:
        Drude/stiffness term from S
        Hall/anomalous Hall term from Omega_z

    second-order response:
        jerk term from J
        nonlinear Hall / BCD term from D

Electric-field convention:
    The script converts E[V/m] to E[meV/nm] using

        E_meV_per_nm = 1e-6 * E_Vm

Thin-film boundary condition:
    E1_t = solve(I + 0.5*sigma1_omega, E1_in)
    source2 = sigma2_omegaomega : E1_t E1_t
    E2_t = -0.5 * solve(I + 0.5*sigma1_2omega, source2)

Environment note:
    Run with a Python environment that has numpy and scipy installed.
    If using conda, activate the correct environment first and run with:

        python transmission_from_response.py ...
"""


from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np

from load_response_mat import load_response_mat


# e^2 / (hbar epsilon0 c), used to convert sheet conductivity in e^2/hbar
# to the dimensionless thin-film boundary-condition convention.
ALPHA_SHEET = 0.09170123610074787


def _normalized_polarization(polarization: tuple[float, float] | np.ndarray) -> np.ndarray:
    pol = np.asarray(polarization, dtype=complex).reshape(2)
    norm = np.linalg.norm(pol)
    if norm == 0:
        raise ValueError("polarization vector must be nonzero")
    return pol / norm


def build_conductivities_from_response(
    S: np.ndarray,
    Omega_z: float,
    D: np.ndarray,
    J: np.ndarray,
    omega_meV: float,
    gamma_meV: float,
) -> dict[str, np.ndarray]:
    """Build dimensionless thin-film conductivities from response tensors.

    Stored response-object dimensions:
      S[gamma,j]      in meV
      Omega_z         dimensionless
      D[i]            in nm
      J[i,j,gamma]    in meV nm
      omega_meV       hbar*omega in meV
      gamma_meV       Gamma in meV

    Returned tensors:
      sigma1_omega[gamma,j]
      sigma1_2omega[gamma,j]
      sigma2_omegaomega[gamma,i,j]

    Dimensions:
      sigma1 is dimensionless.
      sigma2 has units nm/meV and acts on E measured in meV/nm.

    Hall and BCD signs follow the saved response-object convention:
      Omega_z = int f Omega_xy
      D_i = - int Omega_xy partial_i f
      eps = [[0,1],[-1,0]]
    """
    S = np.asarray(S, dtype=complex).reshape(2, 2)
    Omega_z = complex(Omega_z)
    D = np.asarray(D, dtype=complex).reshape(2)
    J = np.asarray(J, dtype=complex).reshape(2, 2, 2)

    eps = np.array([[0, 1], [-1, 0]], dtype=complex)
    z1 = omega_meV + 1j * gamma_meV
    z2 = 2 * omega_meV + 1j * gamma_meV

    sigma1_drude_omega = ALPHA_SHEET * 1j * S / z1
    sigma1_drude_2omega = ALPHA_SHEET * 1j * S / z2
    sigma1_hall = -ALPHA_SHEET * Omega_z * eps

    sigma1_omega = sigma1_drude_omega + sigma1_hall
    sigma1_2omega = sigma1_drude_2omega + sigma1_hall

    # Saved J has order J[i,j,gamma]; the source tensor uses [gamma,i,j].
    sigma2_jerk = ALPHA_SHEET * np.transpose(J, (2, 0, 1)) / (z2 * z1)

    sigma2_bcd = np.zeros((2, 2, 2), dtype=complex)
    for gamma_idx in range(2):
        for i_idx in range(2):
            for j_idx in range(2):
                sigma2_bcd[gamma_idx, i_idx, j_idx] = (
                    ALPHA_SHEET
                    * 1j
                    / z1
                    * 0.5
                    * (
                        eps[gamma_idx, i_idx] * D[j_idx]
                        + eps[gamma_idx, j_idx] * D[i_idx]
                    )
                )

    sigma2_omegaomega = sigma2_jerk + sigma2_bcd

    return {
        "sigma1_drude_omega": sigma1_drude_omega,
        "sigma1_hall": sigma1_hall,
        "sigma1_omega": sigma1_omega,
        "sigma1_drude_2omega": sigma1_drude_2omega,
        "sigma1_2omega": sigma1_2omega,
        "sigma2_jerk": sigma2_jerk,
        "sigma2_bcd": sigma2_bcd,
        "sigma2_omegaomega": sigma2_omegaomega,
    }


def transmitted_harmonics_from_response_tensors(
    S: np.ndarray,
    Omega_z: float,
    D: np.ndarray,
    J: np.ndarray,
    omega_meV: float,
    gamma_meV: float,
    E0_Vm: float,
    polarization: tuple[float, float] | np.ndarray = (1, 0),
) -> dict[str, np.ndarray | complex | float]:
    """Compute transmitted fundamental and second harmonic fields.

    Electric-field conversion:
      e E[V/m] * 1 nm = E * 1e-9 eV = E * 1e-6 meV,
    so E_meV_per_nm = 1e-6 * E_Vm.

    Thin-film boundary condition:
      E1_t = solve(I + 0.5 * sig1_omega, E1_in)
      source2 = einsum("gij,i,j->g", sig2_omegaomega, E1_t, E1_t)
      E2_t = -0.5 * solve(I + 0.5 * sig1_2omega, source2)
    """
    conductivities = build_conductivities_from_response(
        S, Omega_z, D, J, omega_meV, gamma_meV
    )

    pol = _normalized_polarization(polarization)
    identity = np.eye(2, dtype=complex)

    E1_in_Vm = E0_Vm * pol
    E1_in_meVnm = 1e-6 * E1_in_Vm

    E1_t_meVnm = np.linalg.solve(
        identity + 0.5 * conductivities["sigma1_omega"], E1_in_meVnm
    )
    source2_meVnm = np.einsum(
        "gij,i,j->g", conductivities["sigma2_omegaomega"], E1_t_meVnm, E1_t_meVnm
    )
    E2_t_meVnm = -0.5 * np.linalg.solve(
        identity + 0.5 * conductivities["sigma1_2omega"], source2_meVnm
    )

    E1_t_Vm = 1e6 * E1_t_meVnm
    E2_t_Vm = 1e6 * E2_t_meVnm

    E1_norm = np.linalg.norm(E1_t_Vm)
    E2_over_E1_norm = np.linalg.norm(E2_t_Vm) / E1_norm if E1_norm != 0 else np.nan

    return {
        **conductivities,
        "E1_in_Vm": E1_in_Vm,
        "E1_t_Vm": E1_t_Vm,
        "E2_t_Vm": E2_t_Vm,
        "E2_over_E1_norm": E2_over_E1_norm,
    }


def compute_harmonics(
    mat_path: str | Path,
    omega_meV: float,
    mu_meV: float,
    E0_Vm: float,
    polarization: tuple[float, float] | np.ndarray = (1, 0),
    gamma_meV: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load a standardized response object, interpolate at mu, and scatter."""
    response = load_response_mat(mat_path, mu_meV)
    harmonics = transmitted_harmonics_from_response_tensors(
        response["S"],
        response["Omega_z"],
        response["D"],
        response["J"],
        omega_meV=omega_meV,
        gamma_meV=gamma_meV,
        E0_Vm=E0_Vm,
        polarization=polarization,
    )

    result: dict[str, Any] = {
        "material_name": response["material_name"],
        "model_name": response["model_name"],
        "mu_meV": mu_meV,
        "omega_meV": omega_meV,
        "gamma_meV": gamma_meV,
        "S": response["S"],
        "Omega_z": response["Omega_z"],
        "D": response["D"],
        "J": response["J"],
        **harmonics,
    }

    if verbose:
        print(f"Loaded material: {result['material_name']} / {result['model_name']}")
        print(f"Chosen mu: {mu_meV} meV")
        print(f"Chosen omega: {omega_meV} meV")
        print(f"Chosen Gamma: {gamma_meV} meV")
        print(f"E1_in_Vm: {result['E1_in_Vm']}")
        print(f"E1_t_Vm: {result['E1_t_Vm']}")
        print(f"E2_t_Vm: {result['E2_t_Vm']}")
        print(f"E2_over_E1_norm: {result['E2_over_E1_norm']}")

    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute thin-film transmitted harmonics from a response-object .mat file."
    )
    parser.add_argument("response_mat", type=Path)
    parser.add_argument("--omega-meV", type=float, required=True, help="hbar*omega in meV")
    parser.add_argument("--mu-meV", type=float, required=True, help="chemical potential in meV")
    parser.add_argument("--gamma-meV", type=float, default=1.0, help="relaxation scale Gamma in meV")
    parser.add_argument("--E0-Vm", type=float, required=True, help="incident field amplitude in V/m")
    parser.add_argument("--pol", type=float, nargs=2, default=(1.0, 0.0), help="x y polarization vector")
    args = parser.parse_args()

    compute_harmonics(
        args.response_mat,
        omega_meV=args.omega_meV,
        mu_meV=args.mu_meV,
        gamma_meV=args.gamma_meV,
        E0_Vm=args.E0_Vm,
        polarization=np.array(args.pol, dtype=float),
        verbose=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
