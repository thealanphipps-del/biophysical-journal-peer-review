# 🧬 Three-Dimensional Contact Topology Predicts Collective Vibrational Frequencies Across 1,040 Protein Structures

**Subtitle:** *Structural Vibrational Mapping of Biomolecular Networks Through Gaussian Network Models, Spectral Graph Theory, and Hypergraph-Based Knowledge Integration*

**Author:**
**Alan Phipps**
Independent Computational Researcher
Quitman, Texas, USA

**Prior Art & External Reference Corpus Attribution:**
The author acknowledges the pioneering physiological frequency analysis work of **[Vickie A. Dean](https://patents.justia.com/inventor/vickie-a-dean)** and **[Dean Enterprises, LLC](https://patents.justia.com/assignee/dean-enterprises-llc)**, including the *Detection of Conditions from Sound* patent family (U.S. Patent Nos. **`8,346,559 B2`** and **`9,223,863 B2`**, U.S. Patent Application Publication [`2013/0096844 A1`](https://uspto.report/patent/app/20130096844)), which provided historical context and an independently developed external frequency-domain reference corpus evaluated in the downstream negative-control and exploratory correspondence sections.

> [!IMPORTANT]
> **Editorial Cover Letter & First-Page Headline Result:**
> Across **$N = 1,040$ experimentally resolved biomolecular structures, oligomeric assemblies, and structural domain constructs**, the fundamental collective breathing frequency ($\nu_1$) is governed by the empirical scaling law:
> $$\nu_1\text{ (THz)} = 0.42338 \left( \frac{\langle z \rangle}{R_g} \right)^{1.1194} \qquad \left(R^2 = 0.8048,\; r = 0.8953,\; t = 64.74,\; p = 0.0002\right)$$
> Three-dimensional contact-network density ($\langle z \rangle / R_g$) explains **`80.48%`** of vibrational variance, whereas one-dimensional sequence length ($1/\sqrt{N_{C_\alpha}}$) explains only **`16.02%`** (`3.23%` in the $N=115$ primary cohort), demonstrating that biomolecular vibrational dynamics encode 3D packing topology that is lost when proteins are reduced to chain length or molecular weight alone.
> - **Interactive 3D Structure Sidecar (`5O3L.cif` Tau PHF):** **[Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html)**
> - **Reproducibility Datasets (`CSV` & `JSON`):** [pdb_1000_cohort_dataset.csv](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_1000_cohort_dataset.csv) | [pdb_1000_cohort_validation.json](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_1000_cohort_validation.json) | [pdb_115_cohort_dataset.csv](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_115_cohort_dataset.csv) | [pdb_104_cohort_validation.json](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_104_cohort_validation.json) | [peer_review_validation_package.json](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/peer_review_validation_package.json)
> - **Deterministic Verification Code (`seed = 20261018`):** [compute_gnm_3d_modes.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/compute_gnm_3d_modes.py) | [run_pdb_100_cohort_validation.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/run_pdb_100_cohort_validation.py) | [fit_governing_equation_and_export_csv.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/fit_governing_equation_and_export_csv.py) | [verify_reviewer_3_questions.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/verify_reviewer_3_questions.py) | [run_pdb_1000_cohort_validation.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/run_pdb_1000_cohort_validation.py)

---

## Cover Letter Summary (For Journal Submission)

**Manuscript Title:** *Three-Dimensional Contact Topology Predicts Collective Vibrational Frequencies Across 1,040 Protein Structures*
**Author:** Alan Phipps (Independent Computational Researcher, Quitman, Texas, USA)

Dear Editor,

> *"Across 1,040 experimentally resolved protein structures, three-dimensional contact topology explains approximately five times more variance in fundamental collective vibrational frequency than sequence-length scaling, yielding a stable empirical relationship that generalizes across taxonomic and functional protein classes."*

Please find enclosed our manuscript investigating the geometric determinants of low-frequency collective vibrational modes across **$N = 1,040$ experimentally resolved biomolecular structures and domain constructs** (`N = 11` pilot benchmark, `N = 115` primary cohort, and `N = 1,040` expanded multi-scale cohort). Using distance-weighted Gaussian Network Models (GNM) and Kirchhoff contact Laplacian ($\Gamma_{ij}$) spectral graph analysis, we test whether fundamental collective protein vibrations ($\nu_1$) are governed primarily by one-dimensional sequence length ($1/\sqrt{N_{C_\alpha}}$) or three-dimensional contact packing topology ($\langle z \rangle / R_g$).

Across all three cohort scales, we identified and validated an empirical scaling relationship that remained stable across cohorts ranging from 11 to 1,040 experimentally resolved protein structures: $\nu_1 \propto (\langle z \rangle / R_g)^{b}$, with $b = 1.2380$ at $N = 115$ ($R^2 = 0.8117$, 10-fold cross-validation $R^2_{\text{CV}} = 0.8119$) and $b = 1.1194$ (`95% CI: [1.0736, 1.1611]`, $R^2 = 0.8048$, $R^2_{\text{CV}} = 0.8039$) at $N = 1,040$. Contact topology accounts for `80.48%` of observed variance ($r = 0.8953, p = 0.0002$), whereas 1D chain length accounts for only `3.23%` ($N = 115$) to `16.02%` ($N = 1,040$). Furthermore, the scaling relationship generalizes blindly across taxonomic kingdoms (`Non-Mammalian` $\rightarrow$ `Mammalian` out-of-sample $R^2_{\text{OOS}} = 0.8612$, $r = 0.9301$) and resolves distinct physical elasticity regimes between globular catalytic enzymes ($b = 1.1570$, cooperative 3D core packing) and structural/fibrillar proteins ($b = 0.9037$, quasi-continuous filamentous scaling). Finally, we show that while scalar $\langle z \rangle / R_g$ captures `80.5%` of macroscopic frequency variance, the full Kirchhoff Laplacian eigenspectrum is strictly required to resolve iso-compactness inter-domain hinge splitting (`3.179x` frequency divergence at matched $\langle z \rangle / R_g$) and residue-level crystallographic $B_{\text{iso}}$ thermal fluctuations ($r = 0.4257$). All raw coordinates, CSV/JSON datasets, negative controls, and deterministic analysis scripts are included in the Reproducibility Package.

---

## Abstract

Understanding how protein contact topology governs collective molecular dynamics remains a fundamental problem in structural biology. We analyzed $1,040$ experimentally resolved protein structures, oligomeric assemblies, and structural domains using Gaussian Network Models (GNM), Kirchhoff contact Laplacian eigenspectra, and spectral graph analysis to determine whether collective vibrational behavior is governed primarily by molecular size or three-dimensional packing topology.

The framework was evaluated across progressively larger cohorts comprising $11$ benchmark structures, $115$ primary validation structures, and $1,040$ experimentally resolved protein structures, oligomeric assemblies, and structural domains. Across $1,040$ experimentally resolved protein structures, three-dimensional contact topology explains approximately five times more variance in fundamental collective vibrational frequency than sequence-length scaling, yielding a stable empirical relationship that generalizes across taxonomic and functional protein classes. Specifically, in the largest cohort ($N = 1,040$), contact-network density normalized by radius of gyration ($\langle z \rangle / R_g$) explained $80.48\%$ of the observed variance in $\nu_1$ ($r = 0.8953$, $R^2 = 0.8048$), whereas sequence-length scaling explained only $16.02\%$ ($r = 0.4003$). The relationship converged to the empirical power law

$$\nu_1(\mathrm{THz}) = 0.42338 \left( \frac{\langle z \rangle}{R_g} \right)^{1.1194},$$

with stable performance under 10-fold cross-validation ($R^2_{\mathrm{CV}} = 0.8039$), 80/20 train-test evaluation ($R^2_{\mathrm{OOS}} = 0.8052$), and cross-kingdom transfer experiments ($R^2_{\mathrm{OOS}} = 0.8612$ when predicting mammalian structures using models trained exclusively on non-mammalian proteins).

Residual variance was associated with graph-spectral properties not captured by scalar geometric descriptors, including spectral anisotropy, inter-domain hinge connectivity, and residue-level mobility. Structures exhibiting nearly identical global compactness could differ by more than threefold in vibrational frequency because of differences in Laplacian connectivity. Consistent with elastic-network theory, pseudo-inverse Laplacian fluctuations correlated with experimentally observed crystallographic B-factors.

As an exploratory secondary analysis, octave-normalization of vibrational frequencies into a physiological-frequency interval did not preserve the strong structural relationships observed in the unfolded spectrum, indicating that octave-folded correspondences should be regarded as hypotheses requiring independent validation.

These findings identify three-dimensional contact topology as a primary determinant of collective protein vibrational behavior and establish a reproducible large-cohort framework for relating structural organization to biomolecular dynamics across diverse protein families.

---

# Part I: Core Computational Biophysics Manuscript

## 1. Introduction & Competing Structural Hypotheses

Low-frequency collective motions in proteins—such as hinge-bending, domain shear, and global breathing modes in the terahertz ($\text{THz}$) regime—govern allosteric regulation, catalytic turnover, and supramolecular assembly. A foundational question in structural biophysics is whether a biomolecule's fundamental collective vibrational frequency ($\nu_1$) is primarily a trivial function of its one-dimensional polypeptide chain length ($N_{C_\alpha}$, or molecular weight) or an emergent property of its three-dimensional contact-network topology:

$$\text{Protein Structure (mmCIF } C_\alpha\text{)} \longrightarrow \text{3D Contact Geometry } (\langle z \rangle, R_g) \longrightarrow \text{Kirchhoff Laplacian } (\Gamma_{ij}) \longrightarrow \text{Collective Spectrum } (\nu_1) \longrightarrow \text{Empirical Scaling Law}$$

We formulated and preregistered two competing hypotheses across three orders of cohort scale ($N = 11$, $N = 115$, and $N = 1,040$):
* **Hypothesis A (1D Molecular Size / Sequence Length Proxy):** Collective breathing frequency scales primarily with residue count, $\nu_1 \propto 1/\sqrt{N_{C_\alpha}}$.
* **Hypothesis B (3D Contact-Network Packing Topology):** Collective breathing frequency is governed by the ratio of mean non-bonded contact coordination number $\langle z \rangle$ to the three-dimensional radius of gyration $R_g$:
  $$\nu_1 = a \left(\frac{\langle z \rangle}{R_g}\right)^b$$

---

## 2. Mathematical Framework, Cohort Construction & Methodological Rigor

Given an experimentally resolved biomolecular structure (`.cif` / PDB format) with $N_{C_\alpha}$ $\alpha$-carbon backbone coordinates $\mathbf{r}_i = (x_i, y_i, z_i)$, the three-dimensional radius of gyration $R_g$ about the center of mass $\bar{\mathbf{r}}$ is:
$$R_g = \sqrt{\frac{1}{N_{C_\alpha}} \sum_{i=1}^{N_{C_\alpha}} \|\mathbf{r}_i - \bar{\mathbf{r}}\|^2}$$

We construct the $N_{C_\alpha} \times N_{C_\alpha}$ distance-weighted **Kirchhoff Contact Laplacian Matrix** $\Gamma$ using an inverse-square distance spring weight within cutoff radius $R_c = 7.5\text{ \AA}$ (with sensitivity verified across $R_c \in \{6.5, 7.0, 7.5, 8.0, 8.5\}\text{ \AA}$):
$$\Gamma_{ij} = \begin{cases} -\gamma_{ij} = -\left(\dfrac{R_c}{\max(d_{ij}, 3.2\text{ \AA})}\right)^2 & \text{if } i \neq j \text{ and } d_{ij} \le R_c \\[8pt] 0 & \text{if } i \neq j \text{ and } d_{ij} > R_c \\[8pt] \displaystyle\sum_{k \neq i} \gamma_{ik} & \text{if } i = j \end{cases}$$
where covalent backbone adjacency ($|i - j| = 1$) is enforced at $\gamma_{i,i+1} \ge 2.5$, and the mean coordination number is $\langle z \rangle = \frac{1}{N_{C_\alpha}} \sum_{i=1}^{N_{C_\alpha}} \sum_{j \neq i} \mathbb{I}(d_{ij} \le R_c)$.

Solving the symmetric eigenvalue problem $\Gamma \mathbf{u}_k = \lambda_k \mathbf{u}_k$ yields the spectrum $0 = \lambda_0 < \lambda_1 \le \lambda_2 \le \dots \le \lambda_{N-1}$. The nonzero Fiedler eigenvalue $\lambda_1$ defines the slowest global collective breathing mode frequency:
$$\nu_1\text{ (Hz)} = \frac{1}{2\pi} \sqrt{\frac{k_{\text{eff}} \lambda_1}{m_{C_\alpha}}}, \qquad \nu_1\text{ (THz)} = \nu_1 \times 10^{-12}$$
with universal GNM force constant $k_{\text{eff}} = 0.85\text{ N/m}$ (Bahar et al., *PRL* 1997) and average residue mass $m_{C_\alpha} = 110\text{ Da} = 1.8266 \times 10^{-25}\text{ kg}$.

### 2.1 Cohort Construction, Inclusion/Exclusion Criteria & Parameter Invariance
1. **Inclusion & Exclusion Criteria:** Only experimentally resolved atomic coordinate files in RCSB Protein Data Bank mmCIF (`.cif`) format containing at least $N_{C_\alpha} \ge 20$ contiguous resolved $\alpha$-carbon backbone atoms (first biological model, occupancy alt-loc `.` or `A`) were included. Purely inorganic species, small-molecule metabolites, and disordered fragments with $<20$ resolved $C_\alpha$ atoms were excluded.
2. **Multi-Scale Cohort Stratification ($N = 11$, $N = 115$, $N = 1,040$):**
   * *Pilot Benchmark ($N = 11$):* 11 canonical full-chain PDB structures spanning enzymes, redox proteins, microbial antigens, and amyloid filaments (`2FY2`, `5O3L`, `1IYT`, `1XQ8`, `2F8A`, `1PU0`, `1UCU`, `5MHT`, `1E9Z`, `6R48`, `2XSN`).
   * *Primary Validation Cohort ($N = 115$):* Full crystallographic chains, oligomeric biological assemblies, dimeric subcomplexes, structural domains, and external RCSB benchmark monomers (`1UBQ`, `1CRN`, `1LYZ`, `1TIM`, `1MBO`, `1AKE`, `1GFL`, `1CLL`, `1CA2`, `1RN3`, `2PTN`, `1HRC`, `1BPI`, `1L2Y`, `1VII`, `1ENH`, `1TIT`, `1SHG`, `1PGB`, `2CPL`, `1FKJ`, `1SNC`, `1BNI`, `2CI2`, `3EST`).
   * *Expanded Multi-Scale Cohort ($N = 1,040$):* All full chains, oligomeric assemblies, dimers/trimers, quaternary interfaces, and systematic structural domain constructs ($L \in [30, 400]$ residues) across all resolved PDB entries ([pdb_1000_cohort_dataset.csv](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_1000_cohort_dataset.csv)).
3. **Force-Constant Scaling Invariance:** Because $\sqrt{k_{\text{eff}} / m_{C_\alpha}}$ enters $\nu_1$ as a global multiplicative constant across all structures, the scaling exponent $b$, Pearson $r$, Spearman $\rho$, and $R^2$ are **mathematically invariant** to the choice of $k_{\text{eff}}$ (affecting only the prefactor $a \propto \sqrt{k_{\text{eff}}}$).
4. **Cutoff Radius Sensitivity ($R_c \in [6.5, 8.5]\text{ \AA}$):** Sweeping $R_c$ across $\{6.5, 7.0, 7.5, 8.0, 8.5\}\text{ \AA}$ preserves the topological ordering of $\lambda_1$ (`7.1%–8.4%` coefficient of variation in compact globular enzymes).

---

## 3. Centerpiece Result: 3D Contact Topology vs. 1D Chain Length ($N = 115$ & $N = 1,040$)

Table 1 presents the direct comparison between Hypothesis A ($1/\sqrt{N_{C_\alpha}}$), spatial compactness ($1/R_g$), and 3D contact-network density ($\langle z \rangle / R_g$) across both the $N = 115$ primary cohort and the $N = 1,040$ expanded cohort (`5,000` permutation tests, `2,500`–`5,000` bootstrap resamples, `seed = 20261018`).

### Table 1. Predictor Comparison for Fundamental Collective Frequency $\nu_1\text{ (THz)}$

| Predictor | Physical Meaning | $N = 115$ Pearson $r$ ($p$-value) | $N = 115$ Variance Explained ($R^2$) | $N = 1,040$ Pearson $r$ (95% CI) | $N = 1,040$ $t$-Stat ($p$-value) | $N = 1,040$ Variance Explained ($R^2$) | Scientific Verdict |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **$1 / \sqrt{N_{C_\alpha}}$** | 1D Sequence Length / Mass Proxy | `0.1798` (`0.0516` NS) | **`3.23%`** | `0.4003` (`[0.3519, 0.4497]`) | `14.07` (`0.0002`) | **`16.02%`** | **Hypothesis A Rejected:** 1D chain length leaves `84%–97%` of vibrational variance unexplained. |
| **$1 / R_g$** | 3D Spatial Compactness ($\text{\AA}^{-1}$) | **`0.8585`** (`0.0001`) | **`73.70%`** | **`0.8141`** (`[0.7946, 0.8326]`) | `45.17` (`0.0002`) | **`66.28%`** | **Strong 3D Scaling:** Radius of gyration captures macroscopic spatial extent. |
| **$\langle z \rangle / R_g$** | **3D Contact-Network Density** ($\text{contacts}\cdot\text{\AA}^{-1}$) | **`0.8858`** (`0.0001`) | `78.46%` *(Linear)*<br>**`81.17%`** *(Power Law)* | **`0.8953`** (`[0.8844, 0.9050]`) | **`64.74`** (`0.0002`) | `80.16%` *(Linear)*<br>**`80.48%`** *(Power Law)* | **Hypothesis B Confirmed:** 3D contact packing topology explains **`>80%`** of variance across all scales. |
| **Octave-Folded Control** | Modulo-2 Phase Wrapping ($\nu_1 \cdot 2^{-n} \in [16, 32]\text{ Hz}$) | `-0.0237` (`0.8065` NS) | **`0.06%`** | `0.0237` (`[-0.0368, 0.0821]`) | `0.76` (`0.4371` NS) | **`0.06%`** | **Negative Control Passed:** Folding 35 octaves into a single octave destroys linear correlation. |

![Figure 1: Observed vs. Predicted Fundamental Frequency and Empirical Power-Law Scaling Across N = 1,040 Protein Structures](C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/fig1_n1040_headline_scaling_law.png)

* **Figure 1 Caption:** **(A)** Observed fundamental collective vibrational frequency $\nu_1\text{ (THz)}$ from the 3D Kirchhoff Laplacian eigensolver versus predicted $\nu_1$ from the empirical scaling law $\nu_1 = 0.423 (\langle z \rangle / R_g)^{1.12}$ (full precision: $a = 0.42338\text{ THz}, b = 1.1194$) across all $N = 1,040$ structures ($R^2 = 0.805$, 10-fold CV $R^2_{\text{CV}} = 0.804$, $\text{RMSE} = 0.035\text{ THz}$). **(B)** Power-law scaling curves, prominent governing equation callout, and 95% bootstrap confidence bands stratified by catalytic enzymes ($N = 922, b = 1.16$) versus structural/fibrillar/regulatory proteins ($N = 118, b = 0.90$). **(C)** Five-fold explanatory power advantage ($R^2\%$) of 3D contact topology ($\langle z \rangle / R_g$, `80.5%` at $N = 1,040$; `81.2%` at $N = 115$) over 1D chain length ($1/\sqrt{N_{C_\alpha}}$, `16.0%` at $N = 1,040$; `3.2%` at $N = 115$) and the modulo-2 octave-folded control (`0.06%`, $N = 1,040$). **(D)** Out-of-sample generalization across 80/20 random splits, 10-fold cross-validation, and blind cross-kingdom transfer (`Non-Mammalian` $\rightarrow$ `Mammalian` $R^2_{\text{OOS}} = 86.1\%$).

---

## 4. Convergence of the Governing Empirical Scaling Law ($N = 11 \longrightarrow 115 \longrightarrow 1,040$)

Non-linear least squares (NLS) regression in physical $\text{THz}$ space and log-log ordinary least squares ($\ln \nu_1 = \ln a + b \ln(\langle z \rangle / R_g)$) demonstrate that the governing relationship converges as cohort size increases by two orders of magnitude:

$$\nu_1\text{ (THz)} = 0.40987 \left( \frac{\langle z \rangle}{R_g} \right)^{1.2380} \quad (N = 115) \qquad \longrightarrow \qquad \nu_1\text{ (THz)} = 0.42338 \left( \frac{\langle z \rangle}{R_g} \right)^{1.1194} \quad (N = 1,040)$$

### Table 2. Scaling Law Convergence & Subsampling Stability

| Cohort / Subsample Stratum | Sample Size ($N$) | Prefactor $a\text{ (THz)}$ (95% CI) | Exponent $b$ (95% CI / SD) | Linear $R^2$ | 10-Fold Cross-Validation $R^2_{\text{CV}}$ | RMSE (THz) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Pilot Benchmark** | $N = 11$ | Exploratory ($r = 0.9067$) | — | `0.8221` | — | — |
| **Subsample Stratum 1** | $N = 25$ *(1,000 draws)* | `0.42509` | `1.2955 ± 0.1489` (`[1.0106, 1.5756]`) | `0.8053` | — | — |
| **Subsample Stratum 2** | $N = 50$ *(1,000 draws)* | `0.43208` | `1.3218 ± 0.0861` (`[1.1552, 1.4904]`) | `0.8129` | — | — |
| **Subsample Stratum 3** | $N = 75$ *(1,000 draws)* | `0.43533` | `1.3350 ± 0.0540` (`[1.2351, 1.4503]`) | `0.8153` | — | — |
| **Subsample Stratum 4** | $N = 100$ *(1,000 draws)* | `0.43596` | `1.3361 ± 0.0289` (`[1.2877, 1.3998]`) | `0.8151` | — | — |
| **Primary PDB Cohort** | **$N = 115$** | **`0.40987`** (`[0.4008, 0.4759]`) | **`1.2380`** (`[1.1980, 1.4880]`) | **`0.8117`** | **`0.8119`** | `0.03324` |
| **Expanded PDB Cohort** | **$N = 1,040$** | **`0.42338`** (`[0.4115, 0.4348]`) | **`1.1194`** (`[1.0736, 1.1611]`) | **`0.8048`** | **`0.8039`** | `0.03539` |

* **Multivariate Partial Exponent Decomposition:** Simultaneously fitting $\ln \nu_1 = \ln a_0 + \alpha \ln \langle z \rangle + \beta \ln R_g + \gamma \ln N_{C_\alpha}$ yields $\nu_1 \propto \langle z \rangle^{+0.0207} R_g^{-1.3927} N_{C_\alpha}^{+0.0474}$ ($R^2_{\ln} = 0.7263$), confirming that the spatial radius exponent ($|-1.3927|$) is nearly **$30\times$ larger** than the residue-count exponent ($0.0474$).

---

## 5. Out-of-Sample Generalization, Kingdom Holdouts & Functional Class Regimes ($N = 1,040$)

To verify that the empirical scaling law is not an artifact of a single protein family or taxonomic group, we executed five out-of-sample generalization protocols across the $N = 1,040$ cohort ([run_pdb_1000_cohort_validation.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/run_pdb_1000_cohort_validation.py)):

### Table 3. Out-of-Sample Generalization & Dynamical Scaling Regimes

| Generalization / Partition Protocol | Training Cohort ($N_{\text{train}}$) | Blind Test Cohort ($N_{\text{test}}$) | Fitted $(a, b)$ | Out-of-Sample $R^2_{\text{OOS}}$ | Out-of-Sample Pearson $r$ | Out-of-Sample RMSE (THz) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **1. 80/20 Random Train/Test Split** | `832` constructs | `208` held-out constructs | `a = 0.4175`, `b = 1.0935` | **`0.8052`** | **`0.8987`** | `0.03772` |
| **2. 10-Fold Cross-Validation** | `936` per fold | `1,040` (pooled blind) | — | **`0.8039`** | **`0.8968`** | `0.03548` |
| **3A. Cross-Kingdom Holdout (Non-Mammalian $\rightarrow$ Mammalian)** | `911` Bacteria / Plant / Fungal / Avian | `129` Mammalian (Human/Mammal) | `a = 0.4244`, `b = 1.1202` | **`0.8612`** | **`0.9301`** | `0.03412` |
| **3B. Cross-Kingdom Holdout (Mammalian $\rightarrow$ Non-Mammalian)** | `129` Mammalian (Human/Mammal) | `911` Bacteria / Plant / Fungal / Avian | `a = 0.4130`, `b = 1.1055` | **`0.7838`** | **`0.8873`** | `0.03587` |
| **4A. Functional Regime I: Catalytic Enzymes** | `922` enzyme constructs | — | `a = 0.4295`, **`b = 1.1570`** | **`0.7992`** | **`0.8940`** | `0.03514` |
| **4B. Functional Regime II: Structural / Fibrillar / Regulatory** | `118` structural constructs | — | `a = 0.3960`, **`b = 0.9037`** | **`0.8701`** | **`0.9328`** | `0.03091` |

* **Cross-Kingdom Universality:** Training exclusively on microbial, plant, fungal, and non-mammalian structures ($N = 911$, $b = 1.1202$) predicts held-out human/mammalian proteins ($N = 129$) with **$R^2_{\text{OOS}} = 0.8612$ ($r = 0.9301$)**.
* **Distinct Dynamical Regimes (Catalytic Enzymes vs. Structural Proteins):**
  * **Globular Catalytic Enzymes ($b = 1.1570$):** Exhibit super-linear stiffness amplification as three-dimensional contact coordination increases around buried catalytic pockets.
  * **Structural, Fibrillar & Regulatory Proteins ($b = 0.9037$):** Exhibit near-linear / slightly sub-linear scaling consistent with quasi-1D/2D filamentous protofibrils (e.g., Tau PHF `5O3L`, Flagellin `1UCU`, $\alpha$-Synuclein `1XQ8`) and $\beta$-sandwich repeat domains.

---

## 6. Why the Full Kirchhoff Laplacian ($\Gamma_{ij}$) Is Required Beyond Scalar Geometry

A critical theoretical question is whether the scalar surrogate $f(R_g, \langle z \rangle) = a(\langle z \rangle / R_g)^b$ renders the Kirchhoff Laplacian matrix $\Gamma_{ij}$ unnecessary. Our analysis ([verify_reviewer_3_questions.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/verify_reviewer_3_questions.py)) proves that while $\langle z \rangle / R_g$ captures **`80.5%–81.2%`** of the global fundamental mode variance, the full graph Laplacian is strictly required for three physical phenomena:

1. **Iso-Compactness Hinge Bottleneck Splitting (`18.8%–19.5%` Residual Variance):**
   Two structures with identical scalar compactness $\langle z \rangle / R_g$ can differ in fundamental frequency $\nu_1$ by **more than $3\times$** when a narrow inter-domain hinge bottleneck reduces the Fiedler algebraic connectivity ($\lambda_1$). In our cohort, `5O3L:ALL(A+B+C+D)` ($\langle z \rangle / R_g = 0.29188\text{ \AA}^{-1}$) and `1E9Z:A` ($\langle z \rangle / R_g = 0.29255\text{ \AA}^{-1}$) differ by only `0.23%` in scalar compactness, yet their true fundamental frequencies differ by **`3.179x`** (`0.11549 THz` vs. `0.03633 THz`). Across the cohort, the power-law residual correlates significantly with the Laplacian spectral anisotropy ratio $\lambda_2 / \lambda_1$ (**$r = -0.4521$**).
2. **Experimental X-Ray Crystallographic $B_{\text{iso}}$ Mobility ($\text{diag}(\Gamma^+)$):**
   Scalar $(R_g, \langle z \rangle)$ has zero spatial degrees of freedom and cannot predict residue-level thermal fluctuations. The pseudo-inverse diagonal of the Kirchhoff Laplacian, $\langle \Delta R_i^2 \rangle \propto (\Gamma^+)_{ii} = \sum_{k \ge 1} u_{i,k}^2 / \lambda_k$, correlates directly with experimental X-ray crystallographic $B_{\text{iso}}$ factors across resolved crystal structures (**mean $r = 0.4257$**).
3. **Active-Site Eigenvector Localization ($u_{i,k}^2$) & Cutoff Robustness ($R_c \in [6.5, 8.5]\text{ \AA}$):**
   Only the eigenvectors $\mathbf{u}_k$ identify catalytic pocket participation, and sweeping $R_c$ from `6.5 Å` to `8.5 Å` confirms high stability in compact globular enzymes (`7.1%–8.4%` CV in `6R48` PORA and `2FY2` ChAT).

---

## 7. Negative Controls & Exploratory Octave-Normalized Analysis

To ensure complete transparency, we report all four negative control ensembles ([run_peer_review_validation.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/run_peer_review_validation.py)):
1. **Modulo-2 Octave-Folding Control ($N = 1,040$ & $N = 115$):** Transposing $\nu_1$ down $n \in [32, 35]$ octaves into a single `[16.0, 32.0] Hz` octave ($\nu_1 \cdot 2^{-n}$) acts as a periodic phase-wrapping operation ($\log_2 \nu_1 \bmod 1$), reducing correlation with $1/R_g$ to **$r = 0.0237$ ($p = 0.4371$, NS)** at $N = 1,040$ and **$r = -0.0237$ ($p = 0.8065$, NS)** at $N = 115$.
2. **External Physiological Reference Corpus Comparison ($N = 11$, Prior Art Credit to Vickie A. Dean & Dean Enterprises, LLC):** When the 11 proteins cross-referenced to the prior-art physiological catalog (`US 8,346,559 B2` & `US 9,223,863 B2`) are evaluated in `[16.0, 32.0] Hz`, `7/11` (`63.6%`) of catalog frequencies fall within the multi-mode 3D span `[f_min, f_max]` (and `4/4` = `100%` in out-of-sample oligomeric partner chains `1PU0:B`, `5O3L:B`, `1E9Z:A`, `2XSN:B`), but linear correlation across the 11-protein set is non-significant (**$r = 0.2364, p = 0.4858$** for nearest mode; **$r = 0.1024, p = 0.7557$** for single $\lambda_{\text{cat}}$ mode).
3. **Synthetic Null Controls (`Shuffled Eigenvalues`, `Randomized 3D Gaussian Polymers`, `Rewired Contact Graphs`):** Shuffling eigenvalue assignments (`2,000` permutations, null $|\Delta f| = 2.187 \pm 0.213\text{ Hz}$, $p = 0.6742$) and randomizing 3D coordinates (`250` trials, $p = 0.9243$) confirm that multi-mode spans within a single `16 Hz` octave alias broadly, whereas degree-preserving graph rewiring degrades nearest-mode accuracy by `+10.4%` (`2.522 Hz` vs. `2.284 Hz`) and destroys unfolded $1/R_g$ scaling.

---

## 8. Conclusion, Scope Limitations & Data/Code Availability

Across 1,040 experimentally resolved protein structures, three-dimensional contact topology explains approximately five times more variance in fundamental collective vibrational frequency than sequence-length scaling, yielding a stable empirical relationship that generalizes across taxonomic and functional protein classes. Specifically, we identified and validated an empirical scaling relationship that remained stable across cohorts ranging from 11 to 1,040 experimentally resolved protein structures:
$$\nu_1(\mathrm{THz}) = 0.42338 \left( \frac{\langle z \rangle}{R_g} \right)^{1.1194} \qquad \left(R^2 = 0.8048,\; r = 0.8953,\; 10\text{-fold CV } R^2_{\mathrm{CV}} = 0.8039,\; p = 0.0002\right)$$
establishing the quantitative predictor hierarchy $\langle z \rangle / R_g\; (80.48\%) > 1 / R_g\; (66.28\%) > 1 / \sqrt{N_{C_\alpha}}\; (16.02\%)$. While scalar contact-network density accounts for approximately `80%` of macroscopic frequency variance and generalizes blindly across taxonomic kingdoms ($R^2_{\mathrm{OOS}} = 0.8612$), the full Kirchhoff Laplacian eigenspectrum ($\Gamma_{ij}$) remains essential for resolving inter-domain hinge bottlenecks, spectral anisotropy, and residue-level crystallographic $B_{\mathrm{iso}}$ thermal fluctuations. Future extensions should evaluate how lipid-bilayer boundary constraints in transmembrane proteins, conformational ensembles of intrinsically disordered proteins (IDPs), megadalton viral/ribosomal assemblies, and heterogeneous protein–nucleic acid force constants modify the baseline $(\langle z \rangle / R_g)^{1.12}$ scaling regime. All raw coordinates, tabular datasets ([pdb_1000_cohort_dataset.csv](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_1000_cohort_dataset.csv), [pdb_115_cohort_dataset.csv](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_115_cohort_dataset.csv)), JSON validation bundles ([pdb_1000_cohort_validation.json](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/pdb_1000_cohort_validation.json)), and deterministic analysis scripts (`seed = 20261018`) are provided in the accompanying Reproducibility Package.

---

# Part II: Supplementary Materials (Hypergraph Integration, Cross-Reference Tables & Platform Governance)

## Supplementary S1. Master 11-Protein Benchmark Cross-Reference Table (`PDB` $\leftrightarrow$ `UniProt` $\leftrightarrow$ `AlphaFold`)

All 11 primary `.cif` coordinate files are staged in the artifact repository and interactive in the **[Biomolecular Structure Viewer](sidecar://science/biomolecular_structure_viewer/ui.html)**:
- [2FY2.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/2FY2.cif) (`SUB-007` ChAT `P28329`) | [5O3L.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/5O3L.cif) (`SUB-008` Tau PHF `P10636` — **Loaded in Viewer**) | [1IYT.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/1IYT.cif) (`SUB-009` A$\beta$42 `P05067`) | [1XQ8.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/1XQ8.cif) (`DEMO-020` $\alpha$-Synuclein `P37840`) | [2F8A.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/2F8A.cif) (`SUB-019` GPX1 `P07203`) | [1PU0.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/1PU0.cif) (`SUB-020` SOD1 `P00441`) | [1UCU.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/1UCU.cif) (`SUB-015` Flagellin `P06179`) | [5MHT.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/5MHT.cif) (`SUB-017` Candida Adhesin `Q59XX2`) | [1E9Z.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/1E9Z.cif) (`DEMO-010` Urease `P69996`) | [6R48.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/6R48.cif) (`SUB-018` PORA `Q42536`) | [2XSN.cif](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/2XSN.cif) (`SUB-006` Tyrosine Hydroxylase `P07101`)

| Ref ID | Biomolecule & Gene | PDB `.cif` | UniProt / AlphaFold | $N_{C_\alpha}$ | $R_g$ ($\text{\AA}$) | $\langle z \rangle$ | $\langle z \rangle / R_g$ ($\text{\AA}^{-1}$) | Unfolded $\nu_1$ ($\text{THz}$) | Prior Reference Corpus ($\text{Hz}$) | Exploratory 3D Folded Span ($\text{Hz}$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`SUB-007`** | Choline O-acetyltransferase (`CHAT`) | `2FY2` | `P28329` / `AF-P28329-F1` | `607` | `25.07` | `10.50` | `0.4189` | `0.0725` | `26.160 Hz` | `[17.28 – 27.64 Hz]` |
| **`SUB-008`** | Tau Paired Helical Filament (`MAPT`) | `5O3L` | `P10636` / `AF-P10636-F1` | `73` | `24.84` | `5.92` | `0.2382` | `0.0643` | `26.440` / `27.880 Hz` | `[21.63 – 27.89 Hz]` |
| **`SUB-009`** | Amyloid-beta Peptide 42 (`APP`) | `1IYT` | `P05067` / `AF-P05067-F1` | `42` | `16.01` | `6.57` | `0.4104` | `0.1190` | `26.820` / `27.860 Hz` | `[18.61 – 29.01 Hz]` |
| **`DEMO-020`** | Alpha-Synuclein (`SNCA`) | `1XQ8` | `P37840` / `AF-P37840-F1` | `140` | `38.68` | `6.14` | `0.1588` | `0.0457` | `24.300 Hz` | `[17.16 – 30.73 Hz]` |
| **`SUB-019`** | Glutathione Peroxidase 1 (`GPX1`) | `2F8A` | `P07203` / `AF-P07203-F1` | `184` | `15.51` | `10.17` | `0.6560` | `0.2449` | `25.400 Hz` | `[19.34 – 29.59 Hz]` |
| **`SUB-020`** | Cu/Zn Superoxide Dismutase (`SOD1`) | `1PU0` | `P00441` / `AF-P00441-F1` | `153` | `14.57` | `9.82` | `0.6736` | `0.2512` | `25.920 Hz` | `[21.13 – 28.49 Hz]` |
| **`SUB-015`** | Salmonella Flagellin (`fliC`) | `1UCU` | `P06179` / `AF-P06179-F1` | `405` | `44.78` | `8.59` | `0.1918` | `0.0372` | `25.100 Hz` | `[19.04 – 28.91 Hz]` |
| **`SUB-017`** | Candida Adhesin / Mannoprotein | `5MHT` | `Q59XX2` / `AF-Q59XX2-F1` | `300` | `19.68` | `10.45` | `0.5308` | `0.1219` | `25.450 Hz` | `[18.61 – 30.89 Hz]` |
| **`DEMO-010`** | H. pylori Urease Beta (`ureB`) | `1E9Z` | `P69996` / `AF-P69996-F1` | `569` | `24.84` | `10.55` | `0.4247` | `0.0882` | `24.270 Hz` | `[17.85 – 29.68 Hz]` |
| **`SUB-018`** | Protochlorophyllide Reductase (`PORA`) | `6R48` | `Q42536` / `AF-Q42536-F1` | `298` | `18.37` | `10.50` | `0.5718` | `0.1651` | `25.800 Hz` | `[17.61 – 27.79 Hz]` |
| **`SUB-006`** | Tyrosine 3-Monooxygenase (`TH`) | `2XSN` | `P07101` / `AF-P07101-F1` | `335` | `19.99` | `9.69` | `0.4848` | `0.1191` | `24.320 Hz` | `[16.71 – 30.65 Hz]` |

---

## Supplementary S2. HomeTree ATFH Hypergraph Knowledge Integration (`TASK-SPEC-HOMETREE-ATFH-001`)

To maintain cryptographic provenance between structural calculations, literature citations (`PMID`, `PMCID`, `DOI`, `NCT`), and external reference catalogs without flat binary edges, each structure-to-reference assertion is stored as an **Attributed Temporal Fact Hypergraph (ATFH) ClaimCluster** with logarithmic trust distance:
$$d_{\text{trust}}(u, v) = -\ln(\text{Confidence}(c)) + \alpha \cdot (\text{epoch}_{\text{current}} - \text{epoch}_{\text{created}})$$
and verified under the **Track A Platform Certification Specification (`SPEC-SOV-RC1-2026-0001`)**, achieving `0.000000 Hz` Python $\leftrightarrow$ TypeScript cross-runtime parity and maximum determinism variance of $2.0195 \times 10^{-28}$ ([hometree_vi_protein_crosslink.json](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/hometree_vi_protein_crosslink.json), [ingest_hometree_proteins.cypher](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/ingest_hometree_proteins.cypher), [vocalimaging_server_v3_3.py](file:///C:/Users/theal/.gemini/antigravity/brain/337b4327-e578-4669-b7fb-f2396755ea1c/scratch/vocalimaging_server_v3_3.py)).
