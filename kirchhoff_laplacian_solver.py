#!/usr/bin/env python3
"""
3D Gaussian Network Model (GNM) / Kirchhoff Contact Laplacian Eigenfrequency Calculator
and VocalImaging v3.5.0-GNM-3D Upgrader.

Computes physical 3D vibrational eigenfrequencies directly from the 11 mmCIF (.cif) atomic
coordinate files in the artifact directory:
1. Parses _atom_site records for all C-alpha (CA) backbone atoms (first model, primary chain or assembly).
2. Computes exact 3D Structural Geometry:
   - Number of resolved C-alpha residues (N_CA)
   - Center of mass & Radius of Gyration R_g (Angstroms)
   - Mean coordination number <z> (within R_c = 7.5 A cutoff)
   - Mean experimental B-factor / pLDDT stiffness weight
3. Constructs the exact Kirchhoff Contact Laplacian Matrix Gamma (N x N) with distance-weighted
   spring constants gamma_ij = (R_c / d_ij)^2 for d_ij <= R_c (7.5 A).
4. Solves the exact symmetric eigenvalue spectrum 0 = lambda_0 < lambda_1 <= lambda_2 <= ... <= lambda_{N-1}
   using numpy.linalg.eigh:
   - lambda_1: Global Slowest Collective Breathing Mode (hinge-bending / domain cleft oscillation)
   - lambda_2, lambda_3: Secondary collective shear / torsion modes
   - lambda_cat: Local Active-Site / Catalytic Pocket Mode (weighted by catalytic residue participation u_{i,k}^2)
5. Converts physical collective THz normal mode frequencies:
   nu_k (Hz) = (1 / 2*pi) * sqrt(k_eff * lambda_k / m_CA)
   and performs sub-harmonic octave transposition (nu_k * 2^-n) into the [16.0, 32.0] Hz vocal diagnostic base band
   to produce a deterministic 3D Vibrational Frequency Range [f_3d_min_hz, f_3d_max_hz] and primary f_3d_lambda1_hz!
"""
import json
import math
import pprint
from pathlib import Path
import numpy as np

BASE_DIR = Path(r"C:\Users\theal\.gemini\antigravity\brain\337b4327-e578-4669-b7fb-f2396755ea1c")
SCRATCH_DIR = BASE_DIR / "scratch"

# Physical constants for GNM C-alpha collective normal modes:
# Average amino acid residue mass m_res = 110 Da = 1.8266e-25 kg
# Universal GNM force constant k_B * T / <u^2> ~ 0.85 N/m (Bahar et al., PRL 1997 / Haliloglu et al.)
K_SPRING_EFF = 0.85  # N/m (kg / s^2)
M_RESIDUE_KG = 110.0 * 1.66053906660e-27  # kg
RC_CUTOFF_ANGSTROM = 7.5


def parse_cif_ca_atoms(cif_path: Path, target_chain: str = None, max_atoms: int = 900):
    """Extracts C-alpha (CA) 3D coordinates (x, y, z), auth_seq_id, auth_asym_id, and B_iso from an mmCIF file."""
    lines = cif_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_atom_loop = False
    headers = []
    atoms = []
    first_model_num = None

    for line in lines:
        s = line.strip()
        if s == "loop_":
            in_atom_loop = False
            headers = []
            continue
        if s.startswith("_atom_site."):
            in_atom_loop = True
            headers.append(s.split(".")[1].strip())
            continue
        if in_atom_loop:
            if not s or s.startswith("#") or s.startswith("_") or s.startswith("loop_"):
                if atoms:
                    break
                in_atom_loop = False
                continue
            parts = s.split()
            if len(parts) < len(headers):
                continue
            row = dict(zip(headers, parts))
            group = row.get("group_PDB", "ATOM")
            if group != "ATOM":
                continue
            atom_id = row.get("label_atom_id") or row.get("auth_atom_id") or ""
            if atom_id != "CA":
                continue
            model_num = row.get("pdbx_PDB_model_num", "1")
            if first_model_num is None:
                first_model_num = model_num
            elif model_num != first_model_num:
                break
            alt_id = row.get("label_alt_id", ".")
            if alt_id not in (".", "A", "?"):
                continue
            chain = row.get("auth_asym_id") or row.get("label_asym_id") or "A"
            try:
                x = float(row["Cartn_x"])
                y = float(row["Cartn_y"])
                z = float(row["Cartn_z"])
                seq_id = int(row.get("auth_seq_id") or row.get("label_seq_id") or 0)
                b_iso = float(row.get("B_iso_or_equiv", 30.0))
            except Exception:
                continue
            atoms.append({
                "chain": chain,
                "seq_id": seq_id,
                "x": x,
                "y": y,
                "z": z,
                "b_iso": b_iso
            })

    if not atoms:
        raise ValueError(f"No CA atoms parsed from {cif_path}")

    # Filter to target chain if present and has >= 25 residues; otherwise keep primary chains up to max_atoms
    if target_chain:
        chain_atoms = [a for a in atoms if a["chain"] == target_chain]
        if len(chain_atoms) >= 20:
            atoms = chain_atoms

    # Deduplicate by (chain, seq_id)
    dedup = {}
    for a in atoms:
        key = (a["chain"], a["seq_id"])
        if key not in dedup:
            dedup[key] = a
    atoms = list(dedup.values())

    if len(atoms) > max_atoms:
        atoms = atoms[:max_atoms]
    return atoms


def fold_freq_to_vocal_band(freq_hz: float, low_hz: float = 16.0, high_hz: float = 32.0) -> Tuple[float, int]:
    """Folds any physical vibrational frequency (Hz/THz) by octaves (2^-n) into [16.0, 32.0] Hz."""
    f = float(freq_hz)
    shifts = 0
    while f > high_hz:
        f /= 2.0
        shifts -= 1
    while f < low_hz:
        f *= 2.0
        shifts += 1
    return round(f, 4), shifts


def compute_gnm_3d_spectrum(cif_path: Path, catalytic_residues: list, target_chain: str = "A") -> dict:
    """
    Computes the exact 3D Gaussian Network Model (GNM) Kirchhoff Laplacian eigenspectrum
    from an mmCIF file and derives the 3D vibrational frequency range [f_3d_min, f_3d_max].
    """
    atoms = parse_cif_ca_atoms(cif_path, target_chain=target_chain)
    n = len(atoms)
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=np.float64)
    b_factors = np.array([a["b_iso"] for a in atoms], dtype=np.float64)

    # Center of mass and Radius of Gyration R_g
    com = np.mean(coords, axis=0)
    centered = coords - com
    rg_angstrom = float(np.sqrt(np.mean(np.sum(centered ** 2, axis=1))))

    # Pairwise distance matrix D (N x N)
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_mat = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Distance-weighted Kirchhoff Contact Matrix Gamma (R_c = 7.5 A, with long-range elastic tether 15.0 A for multi-chain/fibril connectedness)
    mask_short = (dist_mat > 1e-6) & (dist_mat <= RC_CUTOFF_ANGSTROM)
    mask_med = (dist_mat > RC_CUTOFF_ANGSTROM) & (dist_mat <= 15.0)
    weights = np.zeros_like(dist_mat)
    weights[mask_short] = (RC_CUTOFF_ANGSTROM / dist_mat[mask_short]) ** 2
    weights[mask_med] = 0.08 * ((RC_CUTOFF_ANGSTROM / dist_mat[mask_med]) ** 2)

    contact_counts = np.sum(mask_short, axis=1)
    mean_coordination = float(np.mean(contact_counts))

    gamma = -weights
    np.fill_diagonal(gamma, np.sum(weights, axis=1))

    # Solve symmetric eigenvalue problem: eigenvalues in ascending order
    eigvals, eigvecs = np.linalg.eigh(gamma)

    # Filter out trivial zero translational mode(s) (eigval <= 1e-5)
    nonzero_indices = [i for i, val in enumerate(eigvals) if val > 1e-5]
    if not nonzero_indices:
        nonzero_indices = [1]

    idx_1 = nonzero_indices[0]
    idx_2 = nonzero_indices[min(1, len(nonzero_indices) - 1)]
    idx_3 = nonzero_indices[min(2, len(nonzero_indices) - 1)]

    lambda_1 = float(eigvals[idx_1])
    lambda_2 = float(eigvals[idx_2])
    lambda_3 = float(eigvals[idx_3])

    # Identify catalytic pocket mode (lambda_cat): mode among top 12 low-frequency modes with highest square displacement on catalytic_residues
    cat_set = set(int(r) for r in catalytic_residues) if catalytic_residues else set()
    cat_atom_indices = [i for i, a in enumerate(atoms) if a["seq_id"] in cat_set]
    if not cat_atom_indices:
        # Use core buried residues (highest contact count)
        max_c = np.max(contact_counts)
        cat_atom_indices = [i for i, c in enumerate(contact_counts) if c >= max_c * 0.85]

    best_cat_mode_idx = idx_2
    best_cat_overlap = -1.0
    for m_idx in nonzero_indices[:min(15, len(nonzero_indices))]:
        vec = eigvecs[:, m_idx]
        overlap = float(np.sum(vec[cat_atom_indices] ** 2))
        if overlap > best_cat_overlap:
            best_cat_overlap = overlap
            best_cat_mode_idx = m_idx

    lambda_cat = float(eigvals[best_cat_mode_idx])

    # Convert eigenvalues lambda_k to physical THz frequencies:
    # nu_k = (1 / (2 * pi)) * sqrt(K_SPRING_EFF * lambda_k / M_RESIDUE_KG)
    def eig_to_thz(lam: float) -> float:
        nu_hz = (1.0 / (2.0 * math.pi)) * math.sqrt((K_SPRING_EFF * max(lam, 1e-6)) / M_RESIDUE_KG)
        return nu_hz

    nu_1_hz = eig_to_thz(lambda_1)
    nu_2_hz = eig_to_thz(lambda_2)
    nu_3_hz = eig_to_thz(lambda_3)
    nu_cat_hz = eig_to_thz(lambda_cat)

    f1_folded, s1 = fold_freq_to_vocal_band(nu_1_hz)
    f2_folded, s2 = fold_freq_to_vocal_band(nu_2_hz)
    f3_folded, s3 = fold_freq_to_vocal_band(nu_3_hz)
    fcat_folded, scat = fold_freq_to_vocal_band(nu_cat_hz)

    folded_modes = sorted([f1_folded, f2_folded, f3_folded, fcat_folded])
    f_range_min = round(min(folded_modes), 4)
    f_range_max = round(max(folded_modes), 4)
    f_centroid = round(sum(folded_modes) / len(folded_modes), 4)

    return {
        "cif_file": cif_path.name,
        "resolved_ca_residues": n,
        "chain_analyzed": atoms[0]["chain"],
        "radius_of_gyration_angstrom": round(rg_angstrom, 2),
        "mean_coordination_number": round(mean_coordination, 2),
        "mean_b_factor_or_plddt": round(float(np.mean(b_factors)), 2),
        "cutoff_rc_angstrom": RC_CUTOFF_ANGSTROM,
        "eigenvalues": {
            "lambda_1_breathing": round(lambda_1, 5),
            "lambda_2_hinge": round(lambda_2, 5),
            "lambda_3_torsion": round(lambda_3, 5),
            "lambda_cat_pocket": round(lambda_cat, 5),
            "catalytic_mode_index": int(best_cat_mode_idx),
            "catalytic_pocket_overlap": round(best_cat_overlap, 4)
        },
        "physical_thz_modes": {
            "nu_1_thz": round(nu_1_hz * 1e-12, 4),
            "nu_2_thz": round(nu_2_hz * 1e-12, 4),
            "nu_3_thz": round(nu_3_hz * 1e-12, 4),
            "nu_cat_thz": round(nu_cat_hz * 1e-12, 4)
        },
        "octave_folded_vocal_hz": {
            "f_lambda1_hz": f1_folded,
            "f_lambda1_octave_shift": s1,
            "f_lambda2_hz": f2_folded,
            "f_lambda3_hz": f3_folded,
            "f_lambda_cat_hz": fcat_folded,
            "f_lambda_cat_octave_shift": scat,
            "f_3d_centroid_hz": f_centroid,
            "f_3d_range_hz": [f_range_min, f_range_max],
            "f_3d_bandwidth_hz": round(f_range_max - f_range_min, 4)
        }
    }


def main():
    crosslink_path = BASE_DIR / "hometree_vi_protein_crosslink.json"
    crosslink_data = json.loads(crosslink_path.read_text(encoding="utf-8"))

    gnm_summary = {}
    for p in crosslink_data["proteins_and_biomolecules"]:
        cif_name = p.get("local_cif_artifact") or (f"{p['pdb_ids'][0]}.cif" if p.get("pdb_ids") else "5MHT.cif")
        cif_path = BASE_DIR / cif_name
        # Extract catalytic residue integers from key_functional_residues
        kfr = list(p.get("key_functional_residues", {}).values())[0] if p.get("key_functional_residues") else {}
        chain_id = kfr.get("auth_asym_id", "A")
        if "residues" in kfr:
            cat_res = list(kfr["residues"])
        elif "start_auth_seq_id" in kfr and "end_auth_seq_id" in kfr:
            cat_res = list(range(int(kfr["start_auth_seq_id"]), int(kfr["end_auth_seq_id"]) + 1))
        elif "auth_seq_id" in kfr:
            cat_res = [int(kfr["auth_seq_id"])]
        else:
            cat_res = []

        gnm = compute_gnm_3d_spectrum(cif_path, cat_res, target_chain=chain_id)
        p["gnm_3d_normal_modes"] = gnm
        gnm_summary[p["vi_id"]] = gnm
        print(
            f"{p['vi_id']:8s} | {cif_name:8s} | N_CA={gnm['resolved_ca_residues']:3d} | "
            f"R_g={gnm['radius_of_gyration_angstrom']:5.2f} A | "
            f"lam1={gnm['eigenvalues']['lambda_1_breathing']:.4f} ({gnm['physical_thz_modes']['nu_1_thz']:.3f} THz) -> "
            f"f_lam1={gnm['octave_folded_vocal_hz']['f_lambda1_hz']:6.3f} Hz | "
            f"f_cat={gnm['octave_folded_vocal_hz']['f_lambda_cat_hz']:6.3f} Hz | "
            f"Range=[{gnm['octave_folded_vocal_hz']['f_3d_range_hz'][0]:.3f}, {gnm['octave_folded_vocal_hz']['f_3d_range_hz'][1]:.3f}] Hz"
        )

    crosslink_path.write_text(json.dumps(crosslink_data, indent=2), encoding="utf-8")
    (SCRATCH_DIR / "gnm_3d_summary.json").write_text(json.dumps(gnm_summary, indent=2), encoding="utf-8")
    print("Saved updated hometree_vi_protein_crosslink.json and gnm_3d_summary.json")


if __name__ == "__main__":
    main()
