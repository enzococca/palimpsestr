# Peer Community Journal — Submission Guide

Transfer of the PCI Archaeo–recommended preprint to **Peer Community Journal (PCJ)**.
Diamond open access, **no fees**. PCJ adds the cover page + headers/footers in production.

## Key identifiers (copy-paste ready)

| Field | Value |
|-------|-------|
| Article title | palimpsestR: An R Package for the Identification and Probabilistic Decomposition of Archaeological Palimpsests |
| Preprint id for the portal (ver.2, recommended by PCI) | **10.5281/zenodo.20848849** *(bare DOI, no `https://`)* |
| Latest Zenodo version (v3, with badge) | 10.5281/zenodo.21534452 |
| Zenodo concept DOI (all versions) | 10.5281/zenodo.19969276 |
| PCI recommendation DOI | **10.24072/pci.archaeo.101019** |
| Recommendation page | https://archaeo.peercommunityin.org/articles/rec?id=1019 |
| Recommender | Alessio Palmisano — cite as **Palmisano, 2026** *(confirm on rec page "Cite this recommendation as")* |
| Software / data / code (perennial DOI) | 10.5281/zenodo.19881541 |

## The submission already exists
The portal reported *"Submission with this Preprint id already exists"* — so the submission is
already created. **Do not click Create again.** Go to peercommunityjournal.org → your account →
**"My submissions" / "My articles"**, open it, and **complete it** with the files below.

## What to upload / paste

1. **Formatted manuscript (PDF)** → `palimpsestr-paper-PCJ.pdf`
   PCJ-compliant: starts at Introduction, **no cover page / abstract / badge / line numbers /
   page numbers**; mandatory end sections present; PCJ recommendation sentence in the
   Acknowledgements. **Use this file, NOT the Zenodo `palimpsestr-paper-recommended.pdf`.**

2. **Editable source** → `palimpsestr-paper-PCJ.tex` (generated from the R Markdown source).

3. **BibTeX of the article's references** → `references.bib` (18 refs, this article only).

4. **Article illustration** (png) → `article-illustration.png` (VRPG phase map; we own the rights).

5. **Metadata** (title, authors + ORCID, abstract, keywords, funding, conflict) → paste from
   `02-portal-metadata.md`.

6. **Cover letter** (optional but recommended) → `01-cover-letter.pdf`.

7. **Recommendation PDF** → export it yourself from the browser (Anubis blocks bots):
   open the recommendation page → Print → **Save as PDF** → `pci-recommendation-1019.pdf`.

## Presubmission checklist (PCJ section 5) — status
- [x] Correctly formatted PDF (no cover page/headers/line/page numbers) — `palimpsestr-paper-PCJ.pdf`
- [x] Editable source (.tex) — `palimpsestr-paper-PCJ.tex`
- [x] BibTeX of the article's references — `references.bib`
- [x] Illustration (png/jpeg) we own — `article-illustration.png`
- [x] Mandatory sections present: Acknowledgements · Funding · Conflict of interest disclosure ·
      Data, script, code, and supplementary information availability · References
- [x] PCJ recommendation sentence added to Acknowledgements
- [x] **Perennial DOI for data/scripts/code + a README** describing them — concept DOI
      **10.5281/zenodo.19881541**, which resolves to the latest version. The GitHub→Zenodo
      integration archives every release automatically (v0.24.1 = record 21670989,
      29 Jul 2026), so no manual deposit is needed. The code reference (Cocca 2026) is in
      `references.bib` and cited in the data-availability section.
      **Careful:** `…19881542` is *not* the concept DOI, it is the version DOI of v0.12.0.
      It was cited by mistake in the files sent on 28 Jul and corrected by errata on 29 Jul.
- [ ] Recommendation PDF exported (step 7)

## Production round 1 — reference corrections (27–28 Jul 2026)
The PCJ production office (Elisa Chailler) requested four reference fixes. All done; see
`03-reply-to-production.md` for the point-by-point reply and change log. Reply **by email**
(reply-all to Elisa) attaching `palimpsestr-paper-PCJ.tex`, `references.bib` and
`palimpsestr-paper-PCJ.pdf`. Proofs follow.

## Files in this folder
- `00-submission-guide.md` — this guide
- `01-cover-letter.md` / `01-cover-letter.pdf`
- `02-portal-metadata.md`
- `03-reply-to-production.md` — reply + change log for the production round
- `palimpsestr-paper-PCJ.pdf` — formatted manuscript to submit
- `palimpsestr-paper-PCJ.tex` — editable source
- `references.bib` — BibTeX of the 18 references
- `article-illustration.png` — page illustration
- `pci-recommendation-1019.pdf` — *you export this from the browser (step 7)*

## Notes
- The formatted manuscript is generated from `docs/palimpsestr-paper-pcj.Rmd` (a PCJ variant of
  the paper; the recommended/Zenodo version stays in `docs/palimpsestr-paper.Rmd`).
- Official PCJ LaTeX/docx templates (if a full re-typeset is ever requested): https://osf.io/kmwfv
