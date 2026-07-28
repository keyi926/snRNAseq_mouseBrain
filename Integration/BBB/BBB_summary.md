# BBB Integrity Re-analysis (snRNA-seq only)

**Question:** Does repeated ALNP dosing compromise blood–brain barrier (BBB) integrity?
**Design:** Sample1 = AD model (untreated) · Sample2 = AD model + ALNP (treated) · **n = 1 per group**
**Data:** `Integration/mousebrain_integrated_final.rds` (CCA-integrated, 8,614 nuclei, celltype-annotated)
**Script:** `Integration/snRNAseq_mouseBrain_BBB_integrity.R` → outputs in `Integration/BBB/`

Convention throughout: **"Up" / positive delta = higher in Sample2 (treated)** relative to Sample1 (untreated AD).

---

## Bottom line

Across every BBB-relevant transcriptional program, the ALNP-treated sample shows a barrier that is **maintained or improved**, not compromised. The signal is internally consistent: tight-junction and transporter programs are up, the endothelial-activation program is down, and there is no increase in a leak/transcytosis or immune-infiltration signature. **The snRNA-seq data provide supporting (transcriptional) evidence that repeated ALNP dosing does not disrupt the BBB.**

This is corroborative, not proof — see Caveats.

---

## Results

### 1. Vascular / mural composition (`Vascular_composition_summary.csv`)

| | Endothelial n (frac) | Pericyte n (frac) | Pericyte:EC | Ptprc⁺ frac |
|---|---|---|---|---|
| Sample1 (AD) | 57 (1.32%) | 46 (1.07%) | 0.81 | 2.11% |
| Sample2 (AD+ALNP) | 93 (2.16%) | 45 (1.05%) | 0.48 | 2.16% |

- Endothelial nuclei **not lost** after treatment (captured fraction actually higher).
- Pericyte fraction essentially unchanged (1.07% → 1.05%) → **no pericyte dropout** (pericyte loss is an early sign of BBB breakdown; Nation 2019, Montagne 2020).
- The lower pericyte:EC *ratio* is driven by the higher endothelial count, not pericyte loss — interpret with care.
- **Ptprc⁺ (CD45) fraction unchanged** → no increased peripheral-immune / infiltration proxy.

### 2. Endothelial module scores (`Endothelial_module_scores_summary.csv`, `Endothelial_module_scores.svg`)

| Module | mean S1 | mean S2 | Δ (S2−S1) | cell-level Wilcoxon p |
|---|---|---|---|---|
| **TightJunction** (Cldn5,Ocln,Tjp1,Jam2,Cgn,Esam) | 0.044 | 0.237 | **+0.194** | 0.0039 |
| **Transporter** (Slc2a1,Abcb1a,Abcg2,Slco1a4,Mfsd2a) | 0.884 | 1.082 | **+0.199** | 0.118 |
| **Leak** (Plvap,Cav1) | 0.093 | 0.015 | **−0.078** | 0.392 |
| **EndoActivation** (Vcam1,Icam1,Vwf,Angpt2) | 0.069 | −0.004 | **−0.073** | 0.303 |

Direction: barrier programs ↑, leak & activation programs ↓ — all favouring an intact barrier after treatment.

### 3. Per-gene change in endothelium (`Endothelial_DE_S2_vs_S1_BBBpanel.csv`)

- **Tight junction ↑:** Esam +1.02, Cldn5 +1.00, Ocln +0.82, Tjp1 +0.26, Jam2 +0.14
- **Transporter/efflux ↑:** Mfsd2a +0.48, Slc2a1 (Glut1) +0.32, Slco1a4 +0.25, Abcg2 +0.19, Abcb1a +0.07
- **Endothelial activation ↓:** Angpt2 −1.27, Vcam1 −1.09, Icam1 −0.49, Vwf −0.43
- **Leak / caveolar transcytosis ↓:** Cav1 −0.77; Plvap ≈ 0 in both (no fenestration in either condition)

*(Aqp4/Gja1 appear "up" in the endothelial DE at very low pct — this is ambient/astrocyte-endfoot signal, not endothelial biology; not interpreted.)*

---

## Method support (published precedent)

- **sc/snRNA-seq as a BBB readout in AD:** Yang et al., *Nature* 2022 (human brain vascular atlas) — sc/snRNA detects AD-associated loss of BBB genes in brain vasculature.
- **Mural / endothelial marker definitions:** Vanlandewijck et al., *Nature* 2018 (brain vasculature molecular atlas).
- **Transcytosis / leak markers:** Ben-Zvi 2014 *Nature*, Andreone 2017 *Neuron* (Mfsd2a suppresses caveolar transcytosis; Plvap marks a leaky/fenestrated endothelium).
- **Tight-junction sieve:** Nitta et al., *JCB* 2003 (Cldn5 sets BBB size selectivity).
- **Pericyte loss = early BBB breakdown:** Nation 2019 *Nat Med*; Montagne 2020 *Nature*.
- **Module scoring:** Tirosh et al., *Science* 2016 (AddModuleScore).
- **Statistics caveat:** Squair et al., *Nat Commun* 2021 (pseudoreplication inflates single-cell DE p-values).

---

## Caveats (must be stated when using these results)

1. **n = 1 per group** — no biological-replicate statistics are possible; this is descriptive / hypothesis-supporting, not inferential. Cell-level p-values treat cells as replicates (pseudoreplication) and overstate significance; all BBB-panel `p_val_adj` = 1.
2. **snRNA-seq systematically under-captures endothelium** and does **not** capture infiltrating leukocytes or plasma-protein (IgG/fibrinogen) leakage — the classic functional BBB readouts.
3. **Transcriptome ≠ physical barrier function.** Maintained tight-junction mRNA is supportive but does not prove an intact physical barrier.
4. Effect **direction and consistency** across independent gene programs is the strength here; individual gene magnitudes at these cell numbers are noisy.
