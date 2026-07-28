# Reply to PCJ production office (Elisa Chailler) — revision round 1

**Date of request:** 27 July 2026
**Reply to:** elisa.chailler@inrae.fr (reply-all, keeping PCJ + contact in copy)
**Attachments:** `palimpsestr-paper-PCJ.tex` · `references.bib` · `palimpsestr-paper-PCJ.pdf`

---

## Email text (copy-paste)

Subject: Re: Peer Community Journal — corrected files (palimpsestR)

Dear Elisa,

Thank you very much for the careful check of our reference list. Please find attached the
corrected files: the editable source (`palimpsestr-paper-PCJ.tex`), the bibliography
(`references.bib`), and the recompiled PDF without cover page, headers and footers
(`palimpsestr-paper-PCJ.pdf`).

All four points have been addressed; below is a point-by-point account of the changes.

**1. Royer et al. (2023) — SEAHORS**

You are right that the entry was structurally incomplete. On checking the record, however,
Zenodo 10.5281/zenodo.7957154 is not a software deposit but the *recommended preprint*
(Zenodo resource type: "preprint", version v2), so the "(Version x.y.z) [Software]" form
would not be accurate for it. We have therefore replaced it with the definitive published
version of the same work, which is itself a Peer Community Journal article:

> Royer, A., Discamps, E., Plutniak, S., & Thomas, M. (2023). SEAHORS: Spatial Exploration
> of ArcHaeological Objects in R Shiny. *Peer Community Journal*, 3, e55.
> https://doi.org/10.24072/pcjournal.289

If you would prefer a citation to the software itself rather than to the article, the CRAN
package DOI is https://doi.org/10.32614/cran.package.seahors — we are happy to use that
form instead, or to add it alongside the article; please let us know.

**2. Orphan entry — Cocca (2026)**

Now cited in the text, in the "Data, script, code, and supplementary information
availability" section, immediately after the DOI:

> "… and permanently archived on Zenodo with concept DOI
> https://doi.org/10.5281/zenodo.19881542 (Cocca 2026)."

We have also completed the corresponding `.bib` entry with explicit software/version
information (`version = {0.24.0}`, `note = {R package (Version 0.24.0) [Software]}`), for
consistency with the other software citations.

We re-checked the whole list after these edits: all 18 references are now cited in the
text, and no further orphan entries remain.

**3. Ambroise & Govaert — DOI mismatch**

Your suggestion is correct, and the error was ours: our entry had conflated two different
*Journal of Classification* papers. DOI 10.1007/s003579900012 belongs to Webb (1997),
"Radial basis functions for exploratory data analysis"; the paper we actually cite is
Ambroise & Govaert (1996). The entry now reads:

> Ambroise, C., & Govaert, G. (1996). Constrained clustering and Kohonen Self-Organizing
> Maps. *Journal of Classification*, 13(2), 299–313.
> https://doi.org/10.1007/BF01246104

The in-text citation in the "Methodological framework" section (Neighborhood-EM /
hidden-Markov-random-field stratigraphic term) has been changed accordingly from
"(Ambroise and Govaert 1997)" to "(Ambroise and Govaert 1996)".

**4. Crema & Kobayashi (2020) — "J=omon"**

Corrected to "Jōmon". The stray "=" was an escaping artefact of our LaTeX source (the
macron was written as `J\=omon`, which our converter rendered literally); the reference now
uses the proper Unicode character. The `.bib` entry already carried the correct form
(`J{\=o}mon`) and is unchanged.

**Other changes:** none. Apart from the four corrections above (and the corresponding
change of the BibTeX key `ambroise1997` → `ambroise1996`), the text is identical to the
version you received.

We look forward to receiving the proofs. Please do not hesitate to write if anything above
needs to be handled differently.

With best regards,

Enzo Cocca
on behalf of Enzo Cocca, Roberto Montagnetti and Maurizio Cattani

---

## Change log (internal)

| # | Request | Action | Files |
|---|---------|--------|-------|
| 1 | Royer 2023 structurally incomplete | Replaced Zenodo preprint entry with the published PCJ article (3, e55, `10.24072/pcjournal.289`); verified via Crossref that Zenodo 7957154 is a *preprint* record, not software | `references.bib`, both `.Rmd`, `.tex`, `.pdf` |
| 2 | Cocca 2026 orphan | Cited after the DOI in "Data, script, code, and supplementary information availability"; added `version`/`[Software]` to the bib entry | `references.bib`, both `.Rmd` |
| 3 | Ambroise DOI mismatch | Verified via Crossref: `s003579900012` = Webb 1997. Corrected to 1996, 13(2), 299–313, `10.1007/BF01246104`; in-text year 1997 → 1996; bib key renamed | `references.bib`, both `.Rmd` |
| 4 | "J=omon" typo | `J\=omon` → `Jōmon` (pandoc was rendering the LaTeX macron escape literally) | both `.Rmd` |

Applied to **both** manuscript sources so they do not diverge:
`docs/palimpsestr-paper-pcj.Rmd` (PCJ) and `docs/palimpsestr-paper.Rmd` (Zenodo/recommended).

Post-edit verification: 18/18 references cited, 0 orphans; PDF = 40 pages, starts at
"1 Introduction", no page numbers, no headers/footers.
