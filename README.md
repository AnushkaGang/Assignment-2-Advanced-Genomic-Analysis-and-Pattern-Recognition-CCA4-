# CCA4 – Advanced Genomic Analysis & Pattern Recognition

This repository contains Python code for Assignment 2 (CCA4) Bioinformatics.

## Student Details
**Name:** Anushka Gangwar  
**PRN:** 1032233324  
**Course:** TY BTech CSE Panel A  
**Roll No:** 46

---

## Assignment Structure (Modules for each Question)

| File | Question | Topic |
|------|----------|-------|
| core_dna.py | shared | basic DNA utilities (cleaning / reverse complement) |
| gc_analysis.py | Q1 | GC content analysis, sliding windows, GC skew, comparisons |
| composition.py | Q2 | di/tri nucleotide frequencies, CpG islands, codon usage, classification |
| hamming_distance.py | Q3 | Hamming Distance, distance matrices, ASCII visualization |
| motif_finding.py | Q4 | pattern matching: exact (KMP) + fuzzy (mismatches) + motif conservation |
| orf_reading_frames.py | Q5 | six-frame ORF detection + ORF statistics |
| protein_translation.py | Q6 | DNA→Protein translation, longest ORF, amino acid composition |

Tests:


---

## Running

Run any module individually:

```bash
python gc_analysis.py
python composition.py
python hamming_distance.py
python motif_finding.py
python orf_reading_frames.py
python protein_translation.py

Complexity Notes
Type	Typical Time	Notes
sliding windows	O(n)	single pass
Hamming distance matrix	O(n²) worst case	depends on number of sequences
motif matching (KMP)	O(n + m)	exact search
fuzzy matching	O(n*m)	mismatches allowed
ORFs + translation	O(n)	scanning frames
each module has a menu + takes user input.
