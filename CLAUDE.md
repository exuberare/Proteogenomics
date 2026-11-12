# UAPP: project rules for Claude Code

## What this repo is
Universal Agentic Proteogenomics Platform: Nextflow DSL2 + Python platform for
dark proteome discovery across species. First target: dog (Canis lupus
familiaris, taxon 9615), dataset PXD052932. Next: cat (Felis catus, taxon 9685).

## Where things live
- This repo holds code and configuration only.
- Data lives in ~/uapp_data (= D:\uapp_data): raw/, references/, results/, work/.
- Python: conda env `uapp`. Nextflow: ~/.local/bin/nextflow. Java 21.
- Earlier work: conda env `canine_pg` (WSL) and `darkproteome` (Windows anaconda3). Inspect only.

## Hard rules
1. Inspect before modifying. Explain the plan before editing files.
2. Never delete, move or overwrite anything in ~/uapp_data or /mnt/d. Write new versioned folders.
3. Never change search parameters, FDR thresholds, ORF rules or candidate classes without my
   explicit approval. They live in versioned config files; a change creates a new version.
4. No species, accessions or data paths hardcoded in pipeline code; they belong in config/.
5. Never commit data, FASTA, results or licence keys.
6. Test on tiny data in tests/ first. Ask before any run on real data.
7. Ask before installing software or changing conda environments.
8. Never call a candidate "validated": computational evidence only.
9. Never invent numbers. Report figures only from files or run outputs, with the file path.
10. Ask before git push.

## Conventions
- Nextflow DSL2; one process per modules/local/<name>/main.nf, each with an nf-test test.
- Python package in uapp/, process scripts in bin/, tests with pytest.
- Record each design decision in docs/decisions/NNNN-short-title.md.
