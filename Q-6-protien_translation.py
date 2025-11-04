"""
protein_translation.py — Part C Q6: Translation and amino acid analysis (fresh code).

Implements:
- translate_dna(seq, frame=0) -> AA string (stop='*')
- translate_longest_orf(seq) -> (frame_label, start, end, aa)
- aa_composition(aa_seq) -> freq counts + simple properties

Run:
  python protein_translation.py
"""

from __future__ import annotations
from typing import Dict, Tuple, Optional, List
from core_dna import clean_dna, rev_comp

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

HYDROPHOBIC = set("AILMFWVPG")  # simple set for demo
POLAR = set("STNQYC")
POS = set("KRH")
NEG = set("DE")

def translate_dna(seq: str, frame: int = 0, reverse: bool = False) -> str:
    """
    Translate DNA (ATGC) to AA.
    frame: 0/1/2; reverse=True to use reverse complement strand frames.
    Stops are '*'. No start/stop trimming here.
    """
    s = clean_dna(seq)
    if reverse:
        s = rev_comp(s)
    s = s[frame:]
    aa = []
    for i in range(0, len(s)-2, 3):
        cod = s[i:i+3]
        aa.append(CODON.get(cod, "X"))
    return "".join(aa)

def translate_longest_orf(seq: str) -> Tuple[str, int, int, str]:
    """
    Scan six frames, return the longest ORF translation as (frame_label, start_nt, end_nt, aa_seq).
    Here we consider ORFs that start with 'M' and end at the first stop '*'.
    """
    frames = [("F0",0,False),("F1",1,False),("F2",2,False),("R0",0,True),("R1",1,True),("R2",2,True)]
    best = ("",0,0,"")
    bestlen = -1
    for label, fr, rev in frames:
        aa = translate_dna(seq, frame=fr, reverse=rev)
        i = 0
        while i < len(aa):
            if aa[i] == "M":
                j = i+1
                while j < len(aa) and aa[j] != "*":
                    j += 1
                peptide = aa[i:j]
                if len(peptide) > bestlen:
                    best = (label, i*3, j*3, peptide)
                    bestlen = len(peptide)
                i = j + 1
            else:
                i += 1
    return best

def aa_composition(aa_seq: str) -> Dict[str, float]:
    """Return composition and simple property tallies (fraction)."""
    n = len(aa_seq) or 1
    comp: Dict[str,int] = {}
    for a in aa_seq:
        comp[a] = comp.get(a,0)+1
    out = {k: v/n for k,v in comp.items()}
    out["hydrophobic_frac"] = sum(1 for a in aa_seq if a in HYDROPHOBIC)/n
    out["polar_frac"] = sum(1 for a in aa_seq if a in POLAR)/n
    out["positive_frac"] = sum(1 for a in aa_seq if a in POS)/n
    out["negative_frac"] = sum(1 for a in aa_seq if a in NEG)/n
    return out

# --------------- CLI ---------------
if __name__ == "__main__":
    print("Protein Translation (Q6)")
    s = input("Enter DNA: ")
    print("Translate frames:")
    for lab,(fr,rev) in {"F0":(0,False),"F1":(1,False),"F2":(2,False),"R0":(0,True),"R1":(1,True),"R2":(2,True)}.items():
        aa = translate_dna(s, frame=fr, reverse=rev)
        print(lab, ":", aa[:60] + ("..." if len(aa)>60 else ""))
    lab, st, en, pep = translate_longest_orf(s)
    print("\nLongest ORF:", lab, f"nt:{st}-{en}", f"aa_len:{len(pep)}")
    print("AA composition:", aa_composition(pep))
