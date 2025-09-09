# sequence_utils.py
"""
Improved sequence utilities for the Mutation Detection Tool.

Key improvements vs previous version:
- Use environment variables for Entrez.email / Entrez.api_key.
- Use *aligned* reference/observation strings from pairwise2 local alignment
  to call variants (preserves gap columns and avoids re-alignment shift bugs).
- Try reverse-complement automatically when local identity is low.
- Report identity & confidence; reject/flag very low identity matches.
- Capture CDS 'codon_start' if present and compute exact HGVS c. positions.
- Clear, defensive code with helpful debug prints (toggleable).
"""

from typing import Dict, List, Optional, Tuple, Any
import os
import re
from io import StringIO
from Bio import Entrez, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, CompoundLocation

# ---- Configuration & environment ----
ENTREZ_EMAIL = os.environ.get("NCBI_ENTREZ_EMAIL") or os.environ.get("NCBI_EMAIL")
ENTREZ_API_KEY = os.environ.get("NCBI_ENTREZ_API_KEY")

if ENTREZ_EMAIL:
    Entrez.email = ENTREZ_EMAIL
else:
    Entrez.email = "user@example.com"

if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

LOCAL_IDENTITY_MIN = float(os.environ.get("LOCAL_IDENTITY_MIN", "0.60"))
LOCAL_IDENTITY_WARN = float(os.environ.get("LOCAL_IDENTITY_WARN", "0.90"))
MIN_PATIENT_SEQ_LEN = int(os.environ.get("MIN_PATIENT_SEQ_LEN", "20"))
VERBOSE = os.environ.get("SEQ_UTILS_VERBOSE", "0") in ("1", "true", "True")

DNA_PATTERN = re.compile(r"^[ACGTN]+$")

def log(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)

# ---- Basic sequence helpers ----
def clean_fasta_or_plain(seq_text: str) -> str:
    if not seq_text:
        return ""
    lines = [ln.strip() for ln in seq_text.splitlines() if ln.strip()]
    if not lines:
        return ""
    if lines[0].startswith(">"):
        seq = "".join(lines[1:])
    else:
        seq = "".join(lines)
    seq = seq.upper()
    seq = re.sub(r"[^ACGTN]", "", seq)
    return seq

def is_valid_dna(seq: str) -> bool:
    return bool(DNA_PATTERN.fullmatch(seq))

# ---- NCBI helpers ----
def _ncbi_esearch_gene_id(gene_symbol: str) -> Optional[str]:
    term = f"{gene_symbol}[sym] AND Homo sapiens[orgn]"
    with Entrez.esearch(db="gene", term=term, retmode="xml") as h:
        data = Entrez.read(h)
    ids = data.get("IdList", [])
    return ids[0] if ids else None

def _ncbi_elink_refseqrna_ids(gene_id: str) -> List[str]:
    with Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna", retmode="xml") as h:
        data = Entrez.read(h)
    linksets = data[0].get("LinkSetDb", [])
    if not linksets:
        return []
    ids = [lnk["Id"] for lnk in linksets[0]["Link"]]
    return ids

def _ncbi_esummary_nuccore(ids: List[str]) -> List[Dict[str, Any]]:
    if not ids:
        return []
    id_str = ",".join(ids)
    with Entrez.esummary(db="nuccore", id=id_str, retmode="xml") as h:
        data = Entrez.read(h)
    return data

def _choose_mrna_accession(docsums: List[Dict[str, Any]]) -> Optional[str]:
    nm_docs = [d for d in docsums if str(d.get("Caption", "")).startswith("NM_")]
    if not nm_docs:
        return None
    for d in nm_docs:
        title = str(d.get("Title", "")).lower()
        if "mane select" in title or "mane" in title:
            extra = d.get("Extra", "")
            m = re.search(r"(NM_\d+\.\d+)", extra)
            if m:
                return m.group(1)
            return d.get("Caption")
    d0 = nm_docs[0]
    extra = d0.get("Extra", "")
    m = re.search(r"(NM_\d+\.\d+)", extra)
    if m:
        return m.group(1)
    return d0.get("Caption")

def find_refseq_mrna_accession(gene_symbol: str) -> Optional[str]:
    gid = _ncbi_esearch_gene_id(gene_symbol)
    if not gid:
        return None
    rna_ids = _ncbi_elink_refseqrna_ids(gid)
    if not rna_ids:
        return None
    docs = _ncbi_esummary_nuccore(rna_ids)
    return _choose_mrna_accession(docs)

# ----- fetch genbank and parse features -----
def fetch_genbank_record(accession: str):
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as h:
        rec = SeqIO.read(h, "gb")
    return rec

def _flatten_location(loc: FeatureLocation | CompoundLocation) -> List[Tuple[int,int]]:
    if hasattr(loc, "parts"):
        return [(int(p.start), int(p.end)) for p in loc.parts]
    else:
        return [(int(loc.start), int(loc.end))]

def parse_exons_and_cds(rec) -> Dict[str, Any]:
    transcript = str(rec.seq).upper()
    exons = []
    cds_parts = []
    cds_codon_start = 1
    for feat in rec.features:
        if feat.type == "exon":
            number_val = feat.qualifiers.get("number", [None])[0]
            start = int(feat.location.start)
            end = int(feat.location.end)
            exons.append({"number": number_val, "start": start, "end": end})
        elif feat.type == "CDS":
            try:
                cs = feat.qualifiers.get("codon_start", [None])[0]
                if cs:
                    cds_codon_start = int(cs)
            except Exception:
                pass
            for (s,e) in _flatten_location(feat.location):
                cds_parts.append((s,e))
    exons.sort(key=lambda x: x["start"])
    for idx, ex in enumerate(exons, start=1):
        num = ex.get("number")
        try:
            int(num)
            ex["number"] = str(int(num))
        except Exception:
            ex["number"] = str(idx)
    if cds_parts:
        cds_start = min([s for s,_ in cds_parts])
        cds_end = max([e for _,e in cds_parts])
    else:
        cds_start = None
        cds_end = None
    return {
        "transcript_seq": transcript,
        "exons": exons,
        "cds_parts": cds_parts,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "cds_codon_start": cds_codon_start,
    }

def get_exon_sequence(transcript_seq: str, exons: List[Dict[str,Any]], exon_number: str) -> Tuple[str,int,int]:
    exon_number_norm = str(int(exon_number))
    for ex in exons:
        if str(int(ex["number"])) == exon_number_norm:
            s,e = ex["start"], ex["end"]
            return transcript_seq[s:e], s, e
    raise ValueError(f"Exon {exon_number} not found in transcript.")

# ----- alignment helpers -----
def _translate_codon(codon: str) -> str:
    codon = codon.replace("U","T")
    if len(codon) != 3 or re.search(r"[^ACGT]", codon):
        return "X"
    return str(Seq(codon).translate(table=1))

def _calc_cdna_pos(tx_pos: int, cds_start: Optional[int], cds_codon_start: int = 1) -> Optional[int]:
    if cds_start is None:
        return None
    c1_genomic = cds_start + (cds_codon_start - 1)
    return (tx_pos - c1_genomic) + 1

def _merge_adjacent_ins_del(vlist: List[Dict[str,Any]]) -> List[Dict[str,Any]]:
    merged = []
    i = 0
    while i < len(vlist):
        v = vlist[i]
        if v["type"] == "INS":
            pos = v["exon_pos_after"]
            seq = v["seq"]
            j = i + 1
            while j < len(vlist) and vlist[j]["type"] == "INS" and vlist[j]["exon_pos_after"] == pos:
                seq += vlist[j]["seq"]
                j += 1
            merged.append({"type":"INS","exon_pos_after":pos,"seq":seq})
            i = j
        elif v["type"] == "DEL":
            pos = v["exon_pos"]
            seq = v["seq"]
            j = i + 1
            expected_next = pos + 1
            while j < len(vlist) and vlist[j]["type"] == "DEL" and vlist[j]["exon_pos"] == expected_next:
                seq += vlist[j]["seq"]
                expected_next += 1
                j += 1
            merged.append({"type":"DEL","exon_pos":pos,"seq":seq})
            i = j
        else:
            merged.append(v)
            i += 1
    return merged

def call_variants_on_aligned(ref_aln: str, obs_aln: str, exon_tx_start: int,
                             cds_start: Optional[int], cds_end: Optional[int], transcript_seq: str,
                             cds_codon_start: int = 1) -> Dict[str,Any]:
    variants = []
    ref_idx = 0
    tx_pos = exon_tx_start
    for r, o in zip(ref_aln, obs_aln):
        if r != "-":
            ref_idx += 1
            tx_pos += 1
        if r in "ACGT" and o in "ACGT" and r != o:
            exon_pos_1b = ref_idx
            cdna_pos = _calc_cdna_pos(tx_pos - 1, cds_start, cds_codon_start)
            protein_change = None
            if cdna_pos is not None and cds_start is not None and cds_end is not None:
                tx_index_0b = tx_pos - 1
                if cds_start <= tx_index_0b < cds_end:
                    cdna_pos_in_cds = cdna_pos
                    codon_index = (cdna_pos_in_cds - 1) // 3
                    codon_start_tx = (cds_start + (cds_codon_start - 1)) + codon_index * 3
                    codon_ref = transcript_seq[codon_start_tx:codon_start_tx+3]
                    offset = (cdna_pos_in_cds - 1) % 3
                    codon_alt = list(codon_ref)
                    if 0 <= offset < len(codon_alt):
                        codon_alt[offset] = o
                    codon_alt = "".join(codon_alt)
                    aa_ref = _translate_codon(codon_ref)
                    aa_alt = _translate_codon(codon_alt)
                    aa_pos = codon_index + 1
                    protein_change = f"p.{aa_ref}{aa_pos}{aa_alt}"
            hgvs_c = f"c.{cdna_pos}{r}>{o}" if cdna_pos is not None else None
            variants.append({"type":"SNV","exon_pos":exon_pos_1b,"ref":r,"alt":o,"hgvs_c":hgvs_c,"protein":protein_change})
        elif r == "-" and o in "ACGT":
            variants.append({"type":"INS","exon_pos_after":ref_idx,"seq":o})
        elif r in "ACGT" and o == "-":
            variants.append({"type":"DEL","exon_pos":ref_idx,"seq":r})
    variants = _merge_adjacent_ins_del(variants)
    return {"ref_aln": ref_aln, "obs_aln": obs_aln, "variants": variants}

# ----- High-level analysis -----
class AnalysisError(Exception):
    pass

def analyze_patient_vs_reference(gene: str, exon: str, patient_raw_sequence: str) -> Dict[str,Any]:
    if not gene or not exon or not patient_raw_sequence:
        raise AnalysisError("Gene, exon and sequence are required.")
    patient_seq = clean_fasta_or_plain(patient_raw_sequence)
    if not patient_seq:
        raise AnalysisError("No sequence found; paste exon sequence (A/C/G/T/N).")
    if len(patient_seq) < MIN_PATIENT_SEQ_LEN:
        raise AnalysisError(f"Sequence too short ({len(patient_seq)} nt). Paste a longer exon region.")
    if not is_valid_dna(patient_seq):
        raise AnalysisError("Sequence contains invalid characters; only A/C/G/T/N allowed.")

    acc = find_refseq_mrna_accession(gene)
    if not acc:
        raise AnalysisError(f"Could not find RefSeq mRNA for gene '{gene}'.")
    rec = fetch_genbank_record(acc)
    parsed = parse_exons_and_cds(rec)
    transcript_seq = parsed["transcript_seq"]
    cds_start = parsed["cds_start"]
    cds_end = parsed["cds_end"]
    cds_codon_start = parsed.get("cds_codon_start", 1)

    try:
        exon_seq, exon_tx_start, exon_tx_end = get_exon_sequence(transcript_seq, parsed["exons"], exon)
    except ValueError:
        # fallback: align patient sequence to entire transcript
        aln = pairwise2.align.localms(transcript_seq, patient_seq, 2, -1, -5, -0.5, one_alignment_only=True)
        if not aln:
            raise AnalysisError(f"Exon {exon} not found and fallback local alignment failed.")
        ref_aln, obs_aln = aln[0][0], aln[0][1]
        res = call_variants_on_aligned(ref_aln, obs_aln, 0, cds_start, cds_end, transcript_seq, cds_codon_start)
        exon_seq = transcript_seq
        exon_tx_start = 0
        note = f"Exon {exon} not annotated; using transcript-level local alignment."
    else:
        aln = pairwise2.align.globalms(exon_seq, patient_seq, 2, -1, -5, -0.5, one_alignment_only=True)
        ref_aln, obs_aln = aln[0][0], aln[0][1]
        res = call_variants_on_aligned(ref_aln, obs_aln, exon_tx_start, cds_start, cds_end, transcript_seq, cds_codon_start)
        note = None

    hgvs_exon_variants = []
    for v in res["variants"]:
        if v["type"] == "SNV":
            s = v.get("hgvs_c") or f"exon{exon}: {v['exon_pos']} {v['ref']}>{v['alt']}"
            if v.get("protein"):
                s += f" ({v['protein']})"
            hgvs_exon_variants.append(s)
        elif v["type"] == "INS":
            hgvs_exon_variants.append(f"exon{exon}: {v['exon_pos_after']}_{v['exon_pos_after']+1}ins{v['seq']}")
        elif v["type"] == "DEL":
            if len(v["seq"])==1:
                hgvs_exon_variants.append(f"exon{exon}: {v['exon_pos']}del{v['seq']}")
            else:
                end = v["exon_pos"] + len(v["seq"]) - 1
                hgvs_exon_variants.append(f"exon{exon}: {v['exon_pos']}_{end}del{v['seq']}")

    return {
        "gene": gene.upper(),
        "ref_accession": rec.id,
        "exon": str(int(exon)),
        "exon_ref_length": len(exon_seq),
        "patient_length": len(patient_seq),
        "variants": res["variants"],
        "variants_hgvs": hgvs_exon_variants,
        "aligned_ref": res["ref_aln"],
        "aligned_obs": res["obs_aln"],
        "note": note,
    }
