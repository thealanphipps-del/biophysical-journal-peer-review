# RFC-VI-0004: Dual-Track Evidence Governance, `RULE-VI-27` / `RULE-VI-28` Epistemic Boundaries & 7-Layer Substance-to-Frequency Court of Admissibility (`v4.0.0-court-authority`)

```text
[5DHCP-TKT-0275-DOCMESH-SUBSTRATE]
Context Cube Vector: (X: 3, Y: 2, Face: 2) [Face 2: Execution & Settlement]
Authoritative Base-27 Saevcode: [SAEV-AAABENQAQD0WNNWGPXTQRRIBJFT]
Semantic Focus Drift Vector: ||Δv_f|| = 0.977560 (Cosine Drift: 0.477812) relative to thread root
Lineage Pedigree Link: Root(TKT-8001) → ... → Parent(TKT-0274) → Child(TKT-0275)
Release Tags: v4.0.0-court-authority | v4.0.0-milestone4-final
Pairwise Release Merkle Root (SHA-256): F0D307A9884F863CEC2DCC603B2B9FA67B5D80FC4D4E33181A75F789690188D3
```

---

## 1. Governance Memorandum & Executive Summary

- **To:** Project Stakeholders, Review Board, External Auditors, and Scientific Advisors
- **From:** Governance & Certification Review Function
- **Subject:** Adoption of Dual-Track Evidence Governance and `RULE-VI-27` / `RULE-VI-28` Epistemic Boundary Enforcement
- **Reference:** Release Governance Framework `v4.0.0-court-authority`
- **Status:** Ratified Policy Position

This specification formalizes the permanent governance distinction between **software verification outcomes (`Track 1`)** and **scientific / operational validation outcomes (`Track 2`)** within the Spectral Discovery and Audit Platform (`spectral_discovery_and_audit`) and `SovAreAgentn1` (`n1`).

---

## 2. Authoritative Three-Part & Dual-Track Certification Statements

1. **Engineering Statement (`Track 1`):**
   > *"The implementation behaves according to specification and passed its verification suite."*
2. **Statistical Statement (`Track 1`):**
   > *"The observed result exceeded the configured analytical or null-model thresholds."*
3. **Scientific Statement (`Track 2`):**
   > *"The interpretation remains a hypothesis until supported by independent empirical evidence."*

- **Certified Software Claim:**
  > *"The Spectral Discovery and Audit Platform has successfully completed integrity, provenance, reproducibility, chain-of-custody, adversarial, and regression verification tests and conforms to its documented implementation specification."*
- **Scientific Status Statement:**
  > *"Biological, biochemical, biomarker, or clinical interpretations produced by the platform constitute computational hypotheses and candidate associations. Such interpretations require independent experimental validation and are not established solely through software verification."*
- **Companion Infrastructure Invariant (`RULE-VI-28`):**
  > *"Passing integration tests demonstrates conformance to the implementation specification, but does not by itself establish continuous production operability."*

---

## 3. Five-Layer Epistemic Hierarchy & `RULE-VI-27` Boundary

```text
Software Output
      ↓
Layer 1 Observation               (Measured data)
      ↓
Layer 2 Computational Association (Algorithm execution)
      ↓
Layer 3 Statistical Evidence      (Model and null-model evaluation)
      ↓
Layer 4 Scientific Hypothesis     (Interpretation of computational outputs)
      ↓
────────────────────────────────────────────────────────────
RULE-VI-27 BOUNDARY (No Self-Certification by Software)
────────────────────────────────────────────────────────────
      ↓
Layer 5 Biological Conclusion     (Requires External Empirical Evidence)
```

- **Formal Rule (`RULE-VI-27`):**
  > *"No computational result may be promoted from Layer 4 (Scientific Hypothesis) to Layer 5 (Biological Conclusion) without externally referenced empirical evidence."*

---

## 4. Bound DocMesh Invariants & Cryptographic Checksums

| Parameter | Value | Unit | Bound Source File |
| :--- | :--- | :--- | :--- |
| `UNFOLDED_KIRCHHOFF_R2` | `0.8048` | `R^2` | `engine/scale_erasure.py` |
| `MODULO2_FOLDED_R2` | `0.0006` | `R^2` | `engine/scale_erasure.py` |
| `VARIANCE_COLLAPSE_RATIO` | `1341.33` | `ratio` | `engine/scale_erasure.py` |
| `MEYER_WINDOW_SATURATION_PCT` | `68.39` | `percent` | `engine/chain_of_custody.py` |
| `GATE_A_MIN_SNR` | `1.4` | `ratio` | `engine/chain_of_custody.py` |
| `GATE_A_MAX_SLOPE_DB` | `0.15` | `dB` | `engine/chain_of_custody.py` |
| `RULE_VI_26_REPLAY_TOLERANCE` | `1e-12` | `delta` | `engine/chain_of_custody.py` |
| `MAX_AUTHORIZED_SOFTWARE_LAYER` | `LAYER_4_SCIENTIFIC_HYPOTHESIS` | `enum` | `engine/chain_of_custody.py` |
