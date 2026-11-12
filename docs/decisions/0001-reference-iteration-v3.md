# 0001 - v3 (>= 20 aa) is the reference iteration

Date: 2026-11-12
Decision: The 20 aa lncRNA-ORF iteration (v3) is the reference result.
UAPP's first milestone is to reproduce v3 in Nextflow before any change.

v3 reference values (re-counted from files on 2026-11-12):
- Ensembl ROS_Cfam_1.0 release 115; 7,735 lncRNA transcripts
- 86,084 ORFs >= 20 aa (ATG->stop), 83,062 unique
- Target DB 126,587 (43,525 Ensembl + 83,062 ORFs); 2026-10-02 decoy/contaminant file
- FragPipe per tissue, 27 runs (brain, heart, kidney, liver, lung)
- 1,539 ORF PSMs, 80 peptides, 52 ORFs, 15 final candidates

Not used: v1 (Python, >= 50 aa) and v2 (R, >= 50 aa). The master-prompt figures
(194,430 ORFs; 114,309 targets) match no iteration on disk and are not used.
