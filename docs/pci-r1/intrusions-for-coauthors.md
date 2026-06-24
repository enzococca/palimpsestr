# Poggio Gramignano (VRPG) — find segnalati dal modello: scheda per i co-autori

**A:** Roberto Montagnetti, Maurizio Cattani · **Da:** Enzo Cocca · **Contesto:** revisione PCI Archaeo #1019 (round 1)

## Perché vi chiedo questo

Il referee 1 (G. Bilotti) chiede di **discutere archeologicamente le intrusioni rilevate** nel case study: non basta mostrare che il modello le segnala, bisogna dire se quei reperti sono *plausibili* come materiale residuale / intrusivo / funzionale nel loro contesto, dal punto di vista di chi ha scavato e di chi conosce i materiali.

Il modello (fit K=4, modello categoriale + vincolo di Harris + peso tafonomico + **componente di rumore**) assegna a ogni reperto una **probabilità di essere "fuori contesto"** = posterior della componente uniforme di rumore. Soglia di segnalazione: > 0.5. Sono risultati **10 reperti**, tutti con posterior = **1.0** (il modello non li attribuisce a *nessuna* delle 4 macro-fasi).

⚠️ **Nota metodologica importante (da spiegare nel paper):** poiché a Poggio Gramignano la datazione è *legata alla US* (un solo intervallo di date per unità), la classificazione direzionale risulta `in_context` per tutti: nessun reperto può essere "più vecchio/più giovane" del proprio contesto se condivide la data della US. Quindi la segnalazione qui **non è cronologica**: questi reperti emergono perché la loro US occupa una posizione spaziale/tipologica isolata rispetto ai quattro centroidi di fase (materiale atipico in contesti altrimenti omogenei). La diagnostica direzionale residuale/intrusiva diventa informativa solo con datazione tipologica per-reperto.

## Profilo delle quattro macro-fasi (per riferimento)

| Fase | n reperti | Classe dominante | Intervallo date |
|---|---|---|---|
| 1 | 119 | Reperto Ceramico | 324–425 d.C. |
| 2 | 41 | Reperto Ceramico | 650–350 a.C. (pre-romano) |
| 3 | 135 | Reperto Ceramico | 138–550 d.C. |
| 4 | 320 | Reperto Ceramico | 425–450 d.C. (tardo-antico) |

## I 10 reperti segnalati (raggruppati per US)

"Atipico?" = la classe del reperto ≠ classe dominante della fase a cui è assegnato.
Colonna **Valutazione** = da compilare voi (es. *residuale / intrusivo / funzionale / artefatto di registrazione / coerente*).

| US | Fase | reperti (ID) | Materiale | Date US | taf | Atipico? | Valutazione archeologica (da compilare) |
|---|---|---|---|---|---|---|---|
| **US_204** | 3 | VRPG_0430, 0432, 0433 | Reperto Ceramico | 450–550 d.C. | 0.5 | no | |
| **US_204** | 3 | VRPG_0431 | Materiale Da Costruzione | 450–550 d.C. | 0.5 | **sì** | |
| **US_21** | 3 | VRPG_0435 | Reperto Ceramico | 138–179 d.C. | 0.5 | no | |
| **US_21** | 3 | VRPG_0436 | Materiale Da Costruzione | 138–179 d.C. | 0.5 | **sì** | |
| **US_276** | 3 | VRPG_0489 | Reperto Metallico | 180–324 d.C. | 1.0 | **sì** | |
| **US_302** | 3 | VRPG_0543 | Materiale Da Costruzione | 180–324 d.C. | 0.5 | **sì** | |
| **US_302** | 3 | VRPG_0544 | Reperto Anforaceo | 180–324 d.C. | 0.5 | **sì** | |
| **US_173** | 1 | VRPG_0417 | Materiale Da Costruzione | 324–425 d.C. | 1.0 | **sì** | |

(Tutti: `intrusion_prob` = 1.000 · `direction` = in_context · `intrusion_type` = outlier_in_context. Dati grezzi completi in `intrusions-table.csv`.)

## Domande specifiche per voi

1. **US_204, US_21, US_302** — queste US, con materiale misto (ceramica + materiale da costruzione + anfore) e taf 0.5, sono **riempimenti/livellamenti** o accumuli multi-periodo? Il materiale da costruzione e le anfore in questi contesti sono residuali (riportati da contesti precedenti) o pertinenti?
2. **US_276 (reperto metallico, taf 1.0)** e **US_173 (materiale da costruzione, taf 1.0)** — taf alto (= ben conservato/in posto secondo la scheda), ma il modello li isola. Sono reperti genuinamente singolari per il loro contesto, oppure il taf andrebbe rivisto?
3. C'è qualcuno di questi 10 che attribuireste a **mescolanza funzionale** (come US 304, il focolare nel *dolium*) piuttosto che deposizionale?
4. Qualcuno è invece un **artefatto di registrazione** (due eventi sotto un'unica etichetta US)?

Le vostre risposte confluiranno in un paragrafo nel case study (sezione "Intrusion detection" / "Validation") e nella risposta al referee 1, punto R1.7.
