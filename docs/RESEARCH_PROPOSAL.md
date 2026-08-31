# Redocking LSD: a research proposal

**Working title:** Score-reproducible, pose-unvalidated — what a self-redocking
failure implies for comparative rigid docking of tryptamines at 5-HT2A

**Author:** Fabian Alvarez-Primo · **Project:** [`psilocin-5ht2a`](https://github.com/fabzy4L/psilocin-5ht2a)
· **Related:** [`sert-s438t-escitalopram`](https://github.com/fabzy4L/sert-s438t-escitalopram)
· **Draft date:** 2026-08-31

## Abstract

A rigid AutoDock Vina screen of six serotonergic ligands (LSD, psilocybin,
serotonin, 5-MeO-DMT, psilocin, DMT) against 5-HT2A (PDB 6WGT) produces a
highly reproducible affinity ranking — 5-seed SD ≤0.05 kcal/mol per ligand,
an order of magnitude below Vina's typical ~0.5 kcal/mol noise floor. But
self-redocking 6WGT's own co-crystallized ligand (7LD, confirmed here to be
LSD, not 25-CN-NBOH as commonly mis-cited) into the same box reproduces its
blind-screen affinity almost exactly (−10.07 vs. −10.1 kcal/mol) while
missing its crystal pose by 5.2 Å RMSD — a clear fail against the standard
2.0 Å redocking benchmark. Reproducibility and pose accuracy have come
apart: the ranking is not a search-noise artifact, but it is not evidence
the protocol finds the biologically real pose either. This proposal lays
out the flexible-docking and MD work needed to determine whether the
psilocybin > psilocin ordering — and the ranking generally — reflects real
binding energetics or a scoring-function/rigid-receptor artifact.

## Background

Psilocin (psilocybin's active, dephosphorylated metabolite) and LSD are
both 5-HT2A agonists whose subjective effects can be blunted by SSRI
pretreatment — a clinically observed but mechanistically underspecified
interaction between presynaptic reuptake blockade (SERT) and postsynaptic
agonism (5-HT2A). The sibling project `sert-s438t-escitalopram` found the
same class of gap on the SERT side: rigid docking failed to reproduce a
320-fold experimental Ki loss from the S438T mutation (ΔΔG +0.269 kcal/mol,
within noise), motivating a flexible/MD roadmap there. This project asks
the analogous question on the 5-HT2A side, with a more direct validation
lever available: 6WGT's own bound ligand is one of the six comparators.

## Preliminary data

**Pipeline.** 6WGT (cryo-EM) required a non-trivial prep path: the BRIL
fusion (single-chain, in-frame replacement of ICL3, PDB-renumbered
1001–1106) had to be identified and stripped by walking CA residue numbers
rather than trusting SEQRES; ~80 truncated side chains needed PDBFixer
repair; and meeko's `mk_prepare_receptor` (strict RDKit-template matching)
could not parse the repaired structure at all, failing on an
`AtomValenceException` at the true C-terminus. ADFRsuite's
`prepare_receptor` (geometry-based bond building, no template matching)
succeeded where meeko didn't — a generalizable fallback for any
moderate-resolution cryo-EM GPCR structure meeko chokes on, packaged as
`01c_prep_receptor_adfr.py`.

**Blind-screen ranking, single seed:**

| Ligand | Affinity (kcal/mol) |
|---|---:|
| LSD | −10.1 |
| psilocybin | −7.9 |
| serotonin | −7.1 |
| psilocin | −6.9 |
| 5-MeO-DMT | −6.9 |
| DMT | −6.8 |

**Multi-seed confirmation, n=5 seeds (`04_docking_multiseed.py`):**

| Ligand | mean (kcal/mol) | SD |
|---|---:|---:|
| LSD | −10.045 | 0.021 |
| psilocybin | −7.837 | 0.033 |
| serotonin | −7.083 | 0.019 |
| 5-MeO-DMT | −6.982 | 0.018 |
| psilocin | −6.937 | 0.021 |
| DMT | −6.872 | 0.054 |

The ranking, including psilocybin scoring better than psilocin, is stable
to well under 0.1 kcal/mol across independent search seeds.

**Self-redocking validation (`05_validation_redock.py`).** 7LD's crystal
pose was extracted directly from `data/raw/6WGT.pdb` and redocked into the
same box. Top pose: −10.069 kcal/mol (matches the blind-screen LSD score),
RMSD to the crystal pose: 5.2 Å — a fail against the 2.0 Å convention.
Per-atom deviation was checked across all 24 heavy atoms to rule out a
symmetric-group name-matching artifact (LSD's diethylamide has two
chemically equivalent ethyl arms); deviation was roughly uniform
(2.9–8.1 Å) rather than concentrated in 2–3 atoms, consistent with a
genuinely different bound orientation rather than an RMSD-matching
artifact.

## The problem this proposal addresses

Two results that would normally corroborate each other have split apart.
A reproducible ranking is usually read as evidence a docking protocol is
"working." Here, the same protocol, applied to the one ligand where ground
truth exists, converges on the right affinity magnitude through the wrong
geometry. That means the ranking's reproducibility cannot be used as
indirect evidence for its correctness — the two need to be established
independently, and right now only one is.

## Specific aims

**Aim 1 — Determine whether flexibility recovers the crystal pose.**
Re-run the 7LD self-redock with side-chain flexibility enabled around the
orthosteric pocket residues already flagged as near-box (A:89, 135, 136,
146, 147, 348, 350, 351, 356), and separately against a small ensemble of
receptor conformers (short unbiased MD or normal-mode-perturbed
structures) rather than the single rigid 6WGT coordinate set. Pass
criterion: redocked RMSD <2.0 Å for at least one ensemble member.

**Aim 2 — MD-validate the top-ranked psilocin pose and re-examine the
psilocybin ordering.** Embed the best psilocin pose in a POPC bilayer
(CHARMM-GUI → GROMACS, as already scoped for Phase 2), run 100–200 ns with
GPU offload, and track RMSD/RMSF for pose stability plus the Asp3.32 salt
bridge and TM5 Ser/His contacts. Run MM-GBSA rescoring on trajectory
frames for psilocin and psilocybin to test whether the phosphate-driven
psilocybin advantage survives a flexible, solvated environment — the
specific artifact hypothesis already flagged against the rigid result.

**Aim 3 (stretch) — Signaling-bias proxy.** Compare TM6 outward
displacement in the psilocin trajectory against a Gq-biased agonist
reference structure, as a rough active-state conformation proxy.

## Significance

- **Clinical:** a structurally grounded account of psilocin's 5-HT2A
  engagement is a prerequisite for reasoning about SSRI-psychedelic
  interactions — the question motivating the SERT/5-HT2A pairing.
- **Scientific:** demonstrates, with a same-structure ground-truth case
  rather than an assumed one, that rigid-docking reproducibility is not a
  reliable proxy for pose correctness in this receptor class — relevant
  to any comparative tryptamine/ergoline docking study using GPCR
  cryo-EM structures at this resolution.
- **Methodological:** the ADFR receptor-prep fallback and the paired
  multi-seed/self-redock validation harness are reusable directly (both
  now shared with `sert-s438t-escitalopram`) by anyone hitting meeko's
  template-matching wall on a moderate-resolution structure.

## Limitations

- Single starting receptor conformation (6WGT only); no explicit waters
  or lipids in the current prep.
- Vina's scoring function is empirical and not re-parameterized here for
  ergoline/tryptamine chemistry specifically.
- Six-ligand comparator set is small and biased toward structurally
  related tryptamines/ergolines; no negative control (non-binder).
- Psilocybin is a prodrug — its rigid docking score describes a binding
  event unlikely to be physiologically relevant regardless of pose
  accuracy, which the proposal treats as a feature (a known-artifact case
  to calibrate the flexible/MD pipeline against) rather than noise to
  discard.

## Timeline (single consumer GPU, RTX 5060 8GB)

| Phase | Scope | Est. wall time |
|---|---|---:|
| Aim 1 | Flexible/ensemble redocking validation | 1–2 weeks |
| Aim 2 | MD system prep + 100–200 ns production + MM-GBSA | 3–5 weeks |
| Aim 3 | TM6 bias proxy (stretch) | 1–2 weeks |

## Target venue

Primary: bioRxiv preprint, framed as a methods/case-study note (the
score-vs-pose dissociation plus the reusable ADFR/multi-seed/redock
pipeline travel well as a standalone contribution even before Aim 2
lands). Secondary, once MD data exists: a technical note in *Journal of
Chemical Information and Modeling* or *ACS Chemical Neuroscience*.
