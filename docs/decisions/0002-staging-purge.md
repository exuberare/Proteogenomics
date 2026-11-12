# 0002 - Pipeline may purge its own staging files

Date: 2026-11-12
Reason: C: has 27 GB free and D: 409 GB; datasets are hundreds of GB.
Decision: For each batch the pipeline downloads spectra into
~/uapp_data/staging/<dataset>/<batch>/, searches, extracts results, verifies,
then deletes only that staging folder's mzML/.d/mzBIN/zip files.
Kept: psm/peptide/protein tsv, .pin.gz, workflow, params, logs, PRIDE URL + checksum.
Never deleted: anything outside staging/, existing project folders, references, results.
Claude Code keeps its rm deny rule; only the Nextflow PURGE step deletes.
