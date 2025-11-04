"""
orf_reading_frames.py — Part C Q5: Six-frame ORF detection (fresh code).

Implements:
- six_frames(seq) -> list of frames (3 forward + 3 reverse)
- find_orfs_in_frame(seq, frame_offset) -> list of (start, end, aa_len)
- find_all_orfs(seq, min_len_codon=30) -> dict with frames and ORFs
- simple stats (count, max length)

Run:
  python orf_reading_frames.py
"""

from __future__ import annotations
from typing import Dict, List, Tuple
from core_dna import clean_dna, rev_comp

START = {"ATG"}
STOP = {"TAA","TAG","TGA"}

def six_frames(seq: str) -> List[str]:
    """Return [F0,F1,F2,R0,R1,R2] DNA frames."""
    s = clean_dna(seq)
    rc = rev_comp(s)
    frames = [s[i:] for i in (0,1,2)] + [rc[i:] for i in (0,1,2)]
    return frames

def find_orfs_in_frame(frame_seq: str, frame_offset: int = 0, min_len_codon: int = 30) -> List[Tuple[int,int,int]]:
    """
    Scan a single frame DNA string and return ORFs as (start_nt, end_nt, aa_length).
    Coordinates are relative to the ORIGINAL forward strand for F* frames,
    and to the reverse-complement string for R* frames (reported as indices within the given frame_seq).
    """
    s = frame_seq
    orfs: List[Tuple[int,int,int]] = []
    i = 0
    while i+2 < len(s):
        cod = s[i:i+3]
        if cod in START:
            j = i + 3
            while j+2 < len(s):
                stop_c = s[j:j+3]
                if stop_c in STOP:
                    aa_len = (j - i)//3
                    if aa_len >= min_len_codon:
                        orfs.append((i+frame_offset, j+3+frame_offset, aa_len))
                    i = j + 3
                    break
                j += 3
            else:
                i += 3
        else:
            i += 3
    return orfs

def find_all_orfs(seq: str, min_len_codon: int = 30) -> Dict[str, List[Tuple[int,int,int]]]:
    frames = six_frames(seq)
    labels = ["F0","F1","F2","R0","R1","R2"]
    res: Dict[str, List[Tuple[int,int,int]]] = {}
    for lab, fr in zip(labels, frames):
        res[lab] = find_orfs_in_frame(fr, 0, min_len_codon)
    return res

# --------------- CLI ---------------
if __name__ == "__main__":
    print("ORF detection (Q5)")
    s = input("Enter DNA: ")
    minlen = int(input("Min ORF (aa) [30]: ") or "30")
    orfs = find_all_orfs(s, minlen)
    total = 0
    maxlen = 0
    for k,v in orfs.items():
        print(k, "ORFs:", len(v))
        for (st, en, aa) in v[:10]:
            print(f"  {st}-{en} aa={aa}")
            maxlen = max(maxlen, aa)
        total += len(v)
    print(f"TOTAL ORFs: {total}, LONGEST (aa): {maxlen}")
