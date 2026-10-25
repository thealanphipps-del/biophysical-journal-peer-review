# 📄 *Biophysical Journal* (Cell Press / Biophysical Society) Submission Bundle

**Article Type:** Regular Article — *Computational Biophysics / Protein and Nucleic Acid Dynamics*  
**Manuscript Title:** *Three-Dimensional Contact Topology Predicts Collective Vibrational Frequencies Across 1,040 Protein Structures*  
**Running Title (48 / 50 characters max):** `3D Contact Topology Predicts Protein Vibrations`  
**Author:** **Alan Phipps**¹, *  
¹Independent Computational Biophysics Researcher, Quitman, Texas 75783, USA  
*Correspondence: Alan Phipps (`alan.phipps@pqr.info` | ORCID: `0009-0003-9289-854X`, https://orcid.org/0009-0003-9289-854X)  

**Keywords:** Gaussian Network Model (GNM); Kirchhoff Contact Laplacian; Collective Protein Dynamics; Empirical Scaling Law; Radius of Gyration; Contact-Network Topology

---

## 1. Cover Letter to the Editor-in-Chief (*Biophysical Journal*)

**To:** The Editor-in-Chief and Editorial Board, *Biophysical Journal* (Biophysical Society / Cell Press)  
**Subject:** Submission of Regular Article: *"Three-Dimensional Contact Topology Predicts Collective Vibrational Frequencies Across 1,040 Protein Structures"*

Dear Editor,

> *"Across 1,040 experimentally resolved protein structures, three-dimensional contact topology explains approximately five times more variance in fundamental collective vibrational frequency than sequence-length scaling, yielding a stable empirical relationship that generalizes across taxonomic and functional protein classes."*

Please find enclosed our manuscript, *"Three-Dimensional Contact Topology Predicts Collective Vibrational Frequencies Across 1,040 Protein Structures,"* for consideration as a Regular Article in the *Computational Biophysics* and *Proteins* sections of *Biophysical Journal*.

Since the foundational formulation of the Gaussian Network Model (GNM) by Bahar, Atilgan, and Erman (*Phys. Rev. Lett.* 1997; *Biophys. J.* 2001), elastic network models and Kirchhoff contact Laplacians (**Γᵢⱼ**) have served as standard tools in *Biophysical Journal* for characterizing low-frequency collective protein motions. However, a fundamental quantitative question has remained incompletely resolved across large structural cohorts: *To what extent is a biomolecule's slowest collective breathing frequency (**ν₁**) governed simply by its one-dimensional polypeptide chain length (**1 / √N_Cα**), versus its three-dimensional contact-network packing density (**⟨z⟩ / R_g**)?*

By evaluating distance-weighted Kirchhoff contact Laplacians across progressively larger cohorts (**N = 11** benchmark structures, **N = 115** primary validation structures, and **N = 1,040** experimentally resolved protein structures, oligomeric assemblies, and structural domains), we demonstrate that:

1. **Five-Fold Explanatory Advantage of 3D Contact Topology over Sequence Length:** Across **N = 1,040** structures, three-dimensional contact-network density (**⟨z⟩ / R_g**) explains **80.48%** of the variance in the fundamental collective vibrational mode (**r = 0.8953**, **t = 64.74**, **p = 0.0002**), whereas sequence-length scaling (**1 / √N_Cα**) explains only **16.02%** (**3.23%** in the **N = 115** primary cohort), establishing the strict physical predictor hierarchy **⟨z⟩ / R_g > 1 / R_g > 1 / √N_Cα**.
2. **Convergence of the Empirical Scaling Law:** Across two orders of magnitude in cohort size, the relationship converges onto the stable power law:
   **ν₁ (THz) = 0.42338 · (⟨z⟩ / R_g)^1.1194** &nbsp;&nbsp; *(95% CI for b: [1.0736, 1.1611], 10-fold CV R²_CV = 0.8039)*,  
   which is mathematically invariant in **(b, r, ρ, R²)** to the choice of universal force constant **k_eff** and robust across contact cutoff radii **R_c ∈ [6.5, 8.5] Å**.
3. **Out-of-Sample Cross-Kingdom Transfer & Functional Elasticity Regimes:** Models trained exclusively on non-mammalian structures (**N = 911**, **b = 1.1202**) predict held-out mammalian/human proteins (**N = 129**) blind with **R²_OOS = 0.8612** (**r = 0.9301**). Furthermore, functional stratification reveals distinct physical elasticity regimes between globular catalytic enzymes (**N = 922**, super-linear cooperative core packing **b = 1.1570**) and structural/fibrillar proteins (**N = 118**, near-linear filamentous scaling **b = 0.9037**).
4. **Why the Full Kirchhoff Laplacian (Γᵢⱼ) Remains Necessary:** While scalar **⟨z⟩ / R_g** captures **80.5%** of global fundamental-mode variance, the full Kirchhoff Laplacian eigenspectrum is strictly required to resolve iso-compactness inter-domain hinge bottlenecks (**3.179×** frequency splitting at matched **⟨z⟩ / R_g**, spectral anisotropy correlation **r = −0.4521**) and residue-level X-ray crystallographic **B_iso** thermal fluctuations (**r = 0.4257**).

This manuscript is original, is not under consideration elsewhere, and comes with a cryptographically sealed reproducibility package (`Release ID: SOV-BIOPHYS-N1040-RC1-20261018`, `seed = 20261018`) containing all raw `.cif` coordinates, `CSV`/`JSON` datasets, negative controls, and deterministic Python verification scripts.

Thank you for your time and consideration.

Sincerely,

**Alan Phipps**  
Independent Computational Researcher  
Quitman, Texas, USA

---

## 2. Abstract (*Biophysical Journal* Format — 256 Words)

Understanding how protein contact topology governs collective molecular dynamics remains a fundamental problem in structural biophysics. We analyzed **1,040** experimentally resolved protein structures, oligomeric assemblies, and structural domains using Gaussian Network Models (GNM), Kirchhoff contact Laplacian eigenspectra, and spectral graph analysis to determine whether collective vibrational behavior is governed primarily by molecular size or three-dimensional packing topology. The framework was evaluated across progressively larger cohorts comprising **11** benchmark structures, **115** primary validation structures, and **1,040** experimentally resolved protein structures, oligomeric assemblies, and structural domains. Across **1,040** experimentally resolved protein structures, three-dimensional contact topology explains approximately five times more variance in fundamental collective vibrational frequency than sequence-length scaling, yielding a stable empirical relationship that generalizes across taxonomic and functional protein classes. In the largest cohort (**N = 1,040**), contact-network density normalized by radius of gyration (**⟨z⟩ / R_g**) explained **80.48%** of the observed variance in the fundamental collective mode **ν₁** (**r = 0.8953**, **R² = 0.8048**), whereas sequence-length scaling (**1 / √N_Cα**) explained only **16.02%** (**r = 0.4003**). The relationship converged to the empirical power law **ν₁ (THz) = 0.42338 · (⟨z⟩ / R_g)^1.1194**, with stable performance under 10-fold cross-validation (**R²_CV = 0.8039**), 80/20 train-test evaluation (**R²_OOS = 0.8052**), and cross-kingdom transfer (**R²_OOS = 0.8612**). Residual variance was governed by graph-spectral properties not captured by scalar geometry, including inter-domain hinge bottlenecks (**3.179×** frequency splitting at matched compactness), spectral anisotropy (**λ₂ / λ₁**, **r = −0.4521**), and crystallographic **B_iso** fluctuations (**r = 0.4257**). These findings establish three-dimensional contact topology as a primary determinant of collective protein vibrations across diverse protein families.

---

## 3. Statement of Significance (*Biophysical Journal* Mandatory Section — 118 Words, Limit: ≤ 120 Words)

Low-frequency collective vibrations govern allosteric signaling, enzymatic catalysis, and supramolecular assembly in proteins, yet whether these global breathing modes are dictated primarily by polypeptide chain length or three-dimensional contact packing topology has lacked large-cohort quantitative validation. Across **1,040** experimentally resolved protein structures, oligomeric assemblies, and structural domains, we demonstrate that three-dimensional contact-network density (**⟨z⟩ / R_g**) explains **80.48%** of the variance in fundamental collective vibrational frequency—approximately five times more than sequence-length scaling (**16.02%**). The resulting empirical scaling law, **ν₁ (THz) = 0.42338 · (⟨z⟩ / R_g)^1.1194**, generalizes blindly across taxonomic kingdoms (**R²_OOS = 0.8612**), distinguishes globular enzymes (**b = 1.16**) from fibrillar/structural proteins (**b = 0.90**), and isolates the remaining **19.5%** variance to Laplacian hinge connectivity.

---

## 4. Graphical Abstract & Lead Figure 1 (*Biophysical Journal* 2-Column Layout)

* **Figure 1 Legend:** **(A)** Observed fundamental collective vibrational frequency **ν₁ (THz)** from the 3D Kirchhoff Laplacian eigensolver versus predicted **ν₁** from the empirical scaling law **ν₁ = 0.423 · (⟨z⟩ / R_g)^1.12** (full table precision: **a = 0.42338 THz, b = 1.1194**) across all **N = 1,040** structures (**R² = 0.805**, 10-fold CV **R²_CV = 0.804**, **RMSE = 0.035 THz**). **(B)** Power-law scaling curves, governing equation callout, and 95% bootstrap confidence bands stratified by catalytic enzymes (**N = 922, b = 1.16**) versus structural/fibrillar/regulatory proteins (**N = 118, b = 0.90**). **(C)** Five-fold explanatory power advantage (**R²%**) of 3D contact topology (**⟨z⟩ / R_g**, `80.5%` at **N = 1,040**; `81.2%` at **N = 115**) over 1D chain length (**1 / √N_Cα**, `16.0%` at **N = 1,040**; `3.2%` at **N = 115**) and the modulo-2 octave-folded negative control (`0.06%`, **N = 1,040**). **(D)** Out-of-sample generalization across 80/20 random splits, 10-fold cross-validation, and blind cross-kingdom transfer (`Non-Mammalian → Mammalian` **R²_OOS = 86.1%**).

---

## 5. *Biophysical Journal* Section Cross-Walk & Immutable Submission Archive

The main text, equations, tables, and supplementary materials in `hometree_vocal_imaging_protein_mapping.md` map directly into the required *Biophysical Journal* submission sections:

| *Biophysical Journal* Section | Mapped Manuscript Section in `hometree_vocal_imaging_protein_mapping.md` | Key Contents |
| :--- | :--- | :--- |
| **Title, Running Title, Abstract & Statement of Significance** | Sections 1–3 of this bundle | 48-char running title, 256-word Abstract, 118-word Statement of Significance |
| **Introduction** | `Part I, Section 1` | Competing structural hypotheses: Hypothesis A (**1 / √N_Cα**) vs. Hypothesis B (**⟨z⟩ / R_g**) |
| **Materials and Methods** | `Part I, Section 2 & Section 2.1` | 3D Kirchhoff contact Laplacian (**Γᵢⱼ**, **R_c = 7.5 Å**), Fiedler mode **λ₁**, `.cif` cohort criteria (**N = 11, 115, 1,040**), force-constant invariance, and cutoff sensitivity |
| **Results** | `Part I, Sections 3, 4, 5 & 6` | **Table 1** predictor hierarchy, **Figure 1**, **Table 2** convergence (**N = 11 → 115 → 1,040**), **Table 3** out-of-sample & cross-kingdom generalization (**R²_OOS = 0.8612**), enzyme (**b = 1.1570**) vs. structural (**b = 0.9037**) regimes, and why **Γᵢⱼ** is required beyond scalar geometry (`3.179×` hinge splitting, **B_iso r = 0.4257**) |
| **Discussion & Conclusions** | `Part I, Sections 7 & 8` | Negative controls (`0.06%` octave-folded null, shuffled/Gaussian/rewired nulls), scope limitations (membrane proteins, IDPs, megadalton complexes, nucleoproteins) |
| **Data and Code Availability** | `Part I, Section 8` & `release_manifest_sha256.json` | 26 SHA-256 checksummed files (`seed = 20261018`): `pdb_1000_cohort_dataset.csv`, `pdb_1000_cohort_validation.json`, `pdb_115_cohort_dataset.csv`, and verification scripts |
| **Supporting Material** | `Supplementary S1 & S2` | Representative 11-protein benchmark table (`PDB`, `UniProt`, $N_{C_\alpha}$, $R_g$, $\langle z \rangle$, $\nu_1$) and immutable SHA-256 release manifest (`SOV-BIOPHYS-N1040-RC2-20261019`) |

---

## 6. Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work, the author(s) used Google Gemini and Google Antigravity AI-assisted research and coding tools for Python scientific script verification, mathematical typesetting, and manuscript copy-editing. The author(s) reviewed and edited the output as needed and take full responsibility for the content of the published article.
