# Cohort search (2026-10-22)

Queries and screening behind section 2.4 of `docs/ANALYSIS_PLAN.md`.

| File | Content |
|---|---|
| `ena_studies_2026-10-22.tsv` | ENA portal API, `result=study`, query `(study_title="*periodont*" OR study_description="*periodont*") AND (study_title="*subgingival*" OR study_description="*subgingival*")` |
| `ncbi_bioprojects_2026-10-22.tsv` | NCBI BioProject E-utilities, term `(periodontitis OR periodontal) AND (subgingival OR "dental plaque" OR "periodontal pocket")` |
| `amplicon_runs_2026-10-22.tsv` | Amplicon run summary (runs, library strategy, platform, FASTQ GB) of each project with at least 20 amplicon runs, from the ENA run reports |
| `screening_2026-10-22.tsv` | Projects with at least 20 amplicon runs (from the ENA run reports), plus the SMDI paper's datasets, with a decision and reason each |
| `screen_table.py` | Builds the screening table from the three files above; the decisions are recorded in the script. `python screen_table.py > screening_2026-10-22.tsv` reproduces the committed table exactly |

Decisions were made from the titles, the run reports and the BioSample
attributes (NCBI `efetch db=biosample`, EBI BioSamples API), before any of
the new cohorts was downloaded.
