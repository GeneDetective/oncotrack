# ai_report.py
"""
OpenRouter-only AI exec-summary + local PDF assembly module.

Functions:
 - generate_executive_summary(analysis, patient_info, doctor_notes, desired_output_tokens=900)
     -> str (plain text executive summary, 4 paragraphs)

 - generate_pdf_report(analysis, patient_info, doctor_notes, exec_summary, out_path="genomic_report.pdf")
     -> out_path (writes PDF locally)

 - create_full_report_and_save(analysis, patient_info, doctor_notes, out_path="genomic_report.pdf", use_ai=True, exec_summary=None, **ai_kwargs)
     -> out_path (calls AI via OpenRouter then builds PDF, or uses provided exec_summary)
"""
import os
import re
import time
import tempfile
from typing import Any, Dict, Optional

# requests for OpenRouter
try:
    import requests
    REQUESTS_AVAILABLE = True
except Exception:
    requests = None  # type: ignore
    REQUESTS_AVAILABLE = False

# reportlab for PDF generation
try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    REPORTLAB_AVAILABLE = True
except Exception:
    REPORTLAB_AVAILABLE = False

# -------------------------
# OpenRouter settings
# -------------------------
OPENROUTER_API_KEY = os.environ.get("OPENROUTER_API_KEY", None)
OPENROUTER_BASE = os.environ.get("OPENROUTER_BASE", "https://openrouter.ai/api/v1")
# default model name (OpenRouter accepts provider suffixes like ":free")
OPENROUTER_MODEL = os.environ.get("OPENROUTER_MODEL", "openai/gpt-oss-20b:free")

OPENROUTER_AVAILABLE = bool(OPENROUTER_API_KEY and REQUESTS_AVAILABLE)

# small regex helpers
_CONTINUE_TOKEN_RE = re.compile(r"\bCONTINUE_TOKEN\s*:\s*(\S+)", re.IGNORECASE)
_PART_RE = re.compile(r"\bPART\s*(\d+)\s*[/\\]\s*(\d+)\b", re.IGNORECASE)

# ----------------------
# Token estimate helpers
# ----------------------
def estimate_tokens_from_text(text: str) -> int:
    if not text:
        return 0
    # Conservative: 1 token ≈ 3.5 characters
    return max(1, int(len(text) / 3.5))


def choose_max_tokens_for_output(desired_output_tokens: int, prompt_text: str, model_context_limit: int = 8192) -> int:
    prompt_tokens = estimate_tokens_from_text(prompt_text)
    margin = int(prompt_tokens * 0.2) + 50
    available = max(0, model_context_limit - prompt_tokens - margin)
    return min(desired_output_tokens, max(16, available))


# ----------------------
# System prompt & builder
# ----------------------
_SYSTEM_TEXT = (
    "You are a professional clinical genomics reporting assistant for oncologists and molecular pathologists.\n"
    "IMPORTANT (follow exactly): Output ONLY a single plain-text Executive Summary (no markdown, no headings, no tables, no JSON). "
    "Return exactly FOUR short paragraphs separated by one blank line. "
    "Paragraph1: Patient profile (1-2 sentences). Paragraph2: Mutation biology (1-3 sentences). "
    "Paragraph3: Therapeutic implications (2-4 short sentences). Paragraph4: Conclusion & immediate next tests (1-2 sentences). "
    "Do not invent IDs or fabricate trial numbers. Use neutral wording (likely, may, suggests). Keep within the token budget provided by the client."
)


def build_exec_prompt(analysis: Dict, patient_info: Dict, doctor_notes: Optional[str]) -> str:
    gene = (analysis.get("gene") or "N/A").strip()
    exon = analysis.get("exon") or "N/A"
    ref = (analysis.get("ref_accession") or "N/A").strip()
    variants = analysis.get("variants_hgvs") or []
    variants_text = "; ".join(variants) if variants else "No variant detected"
    alignment_ref = (analysis.get("aligned_ref") or "").strip().replace("\n", " ")
    alignment_obs = (analysis.get("aligned_obs") or "").strip().replace("\n", " ")
    alignment_excerpt = ""
    if alignment_ref or alignment_obs:
        a_ref = (alignment_ref[:240] + "...") if len(alignment_ref) > 240 else alignment_ref
        a_obs = (alignment_obs[:240] + "...") if len(alignment_obs) > 240 else alignment_obs
        alignment_excerpt = f"Alignment excerpt (REF / PAT):\nREF: {a_ref}\nPAT: {a_obs}\n"

    user_block = (
        f"PATIENT DATA:\n"
        f"- Age: {patient_info.get('age', 'N/A')}\n"
        f"- Sex: {patient_info.get('sex', 'N/A')}\n"
        f"- Cancer type: {patient_info.get('cancer_type', 'N/A')}\n"
        f"- Doctor notes: {doctor_notes or 'None'}\n\n"
        f"VARIANT CALL:\n"
        f"- Gene: {gene}\n"
        f"- RefSeq mRNA: {ref}\n"
        f"- Exon: {exon}\n"
        f"- Variant(s): {variants_text}\n\n"
        f"{alignment_excerpt}"
        "TASK: Generate ONLY the Executive Summary as plain text (4 paragraphs). Keep it clinician-focused and concise.\n"
    )
    return user_block


# ----------------------
# OpenRouter call helper (matches test script behaviour)
# ----------------------
def _call_openrouter_chat(system_msg: str, user_msg: str, model: str, max_tokens: int, temperature: float = 0.0, timeout: int = 60) -> Dict[str, Any]:
    """
    Call OpenRouter chat completions endpoint via requests.
    Returns parsed JSON (dict) on success; raises RuntimeError with status/body on non-200.
    """
    if not OPENROUTER_AVAILABLE:
        raise RuntimeError("OpenRouter not available: set OPENROUTER_API_KEY and install 'requests'.")

    url = OPENROUTER_BASE.rstrip("/") + "/chat/completions"
    headers = {
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "Content-Type": "application/json",
    }
    # optional extra headers for OpenRouter leaderboard (safe to omit)
    referer = os.environ.get("OPENROUTER_HTTP_REFERER")
    xtitle = os.environ.get("OPENROUTER_X_TITLE")
    if referer:
        headers["HTTP-Referer"] = referer
    if xtitle:
        headers["X-Title"] = xtitle

    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": system_msg},
            {"role": "user", "content": user_msg}
        ],
        "max_tokens": int(max_tokens),
        "temperature": float(temperature),
    }

    resp = requests.post(url, headers=headers, json=payload, timeout=timeout)
    # return full info for inspection / continuation logic
    if resp.status_code != 200:
        # include small excerpt of body for debug
        body = resp.text
        raise RuntimeError(f"OpenRouter API error {resp.status_code}: {body}")
    # If response content-type isn't JSON something went wrong, but try to parse
    return resp.json()


# ----------------------
# Extraction helper (OpenRouter/OpenAI-like response)
# ----------------------
def _extract_text(resp: Any) -> str:
    """
    Robustly extract assistant text from OpenRouter/OpenAI-style response JSON.
    """
    try:
        if isinstance(resp, dict):
            choices = resp.get("choices") or []
            if choices:
                ch = choices[0]
                # OpenRouter style: ch may contain 'message' with 'content' string or dict
                if isinstance(ch, dict):
                    msg = ch.get("message") or ch.get("delta") or {}
                    if isinstance(msg, dict):
                        content = msg.get("content")
                        if isinstance(content, str) and content.strip():
                            return content.strip()
                        if isinstance(content, dict):
                            # try parts/text inside content
                            for k in ("text", "content", "parts", "body"):
                                if k in content:
                                    v = content[k]
                                    if isinstance(v, str) and v.strip():
                                        return v.strip()
                                    if isinstance(v, list) and v:
                                        return " ".join([str(x) for x in v if isinstance(x, str)]).strip()
                # fallback keys
                for k in ("text", "generated_text", "output_text", "reasoning_content"):
                    if isinstance(ch.get(k), str) and ch.get(k).strip():
                        return ch.get(k).strip()
            # top-level fallback strings
            for v in resp.values():
                if isinstance(v, str) and v.strip():
                    return v.strip()
        if isinstance(resp, str):
            return resp.strip()
    except Exception:
        pass
    return ""


# ----------------------
# Public: generate_executive_summary (OpenRouter only) with continuation handling
# ----------------------
def generate_executive_summary(analysis: Dict, patient_info: Dict, doctor_notes: Optional[str] = None,
                               desired_output_tokens: int = 900, model_context_limit: int = 8192,
                               model_id: Optional[str] = None) -> str:
    """
    Generate an executive summary calling OpenRouter.
    model_id (optional) overrides OPENROUTER_MODEL for this call.

    If the first reply is truncated (finish_reason != "stop"), attempt up to 2 continuation calls
    that include the assistant's previous content in the message history so the model can continue.
    """
    if not OPENROUTER_AVAILABLE:
        raise RuntimeError("OpenRouter API key not configured or requests not installed. Set OPENROUTER_API_KEY and pip install requests.")

    model = model_id or OPENROUTER_MODEL
    user_block = build_exec_prompt(analysis, patient_info, doctor_notes)
    # decide how many tokens to request (best-effort)
    max_tokens = choose_max_tokens_for_output(desired_output_tokens, _SYSTEM_TEXT + "\n" + user_block, model_context_limit)
    if max_tokens < 64:
        raise RuntimeError("Not enough context available for generation; shorten inputs.")

    try:
        # first call
        raw = _call_openrouter_chat(_SYSTEM_TEXT, user_block, model=model, max_tokens=max_tokens, temperature=0.0)
        # extract text and check finish_reason
        text = _extract_text(raw)
        finish_reason = None
        try:
            choices = raw.get("choices") or []
            if choices and isinstance(choices, list) and isinstance(choices[0], dict):
                finish_reason = choices[0].get("finish_reason") or choices[0].get("native_finish_reason")
        except Exception:
            finish_reason = None

        # if truncated, attempt to continue (up to 2 times)
        attempts = 0
        combined = text or ""
        # Use the assistant's previous content as context for continuation
        while attempts < 2 and finish_reason and finish_reason.lower() in ("length", "max_tokens", "length; stop", "incomplete"):
            attempts += 1
            # Build messages that include the previous assistant text so model can continue
            # Note: OpenRouter expects a full messages list for each call
            prev_assistant = combined
            continue_user_msg = "Please continue the previous assistant response from where it left off. Keep the same clinical tone and complete the 4-paragraph executive summary."
            # For continuation, call endpoint with messages containing system, user, assistant(prev), user(continue)
            url = OPENROUTER_BASE.rstrip("/") + "/chat/completions"
            headers = {
                "Authorization": f"Bearer {OPENROUTER_API_KEY}",
                "Content-Type": "application/json",
            }
            referer = os.environ.get("OPENROUTER_HTTP_REFERER")
            xtitle = os.environ.get("OPENROUTER_X_TITLE")
            if referer:
                headers["HTTP-Referer"] = referer
            if xtitle:
                headers["X-Title"] = xtitle

            payload = {
                "model": model,
                "messages": [
                    {"role": "system", "content": _SYSTEM_TEXT},
                    {"role": "user", "content": user_block},
                    {"role": "assistant", "content": prev_assistant},
                    {"role": "user", "content": continue_user_msg}
                ],
                # try to request additional tokens -- increase on continuation
                "max_tokens": int(min(max_tokens * 1.5, 2000)),
                "temperature": 0.0,
            }
            resp = requests.post(url, headers=headers, json=payload, timeout=60)
            if resp.status_code != 200:
                # stop trying on provider error; include body
                raise RuntimeError(f"OpenRouter API error during continuation {resp.status_code}: {resp.text}")
            raw2 = resp.json()
            cont_text = _extract_text(raw2)
            # append new text (naive concatenation)
            if cont_text:
                # avoid repeating identical prefix
                if cont_text.startswith(prev_assistant.strip()):
                    combined = cont_text
                else:
                    combined = (combined + "\n\n" + cont_text).strip()
            # update finish_reason for possible further continuation
            try:
                choices2 = raw2.get("choices") or []
                if choices2 and isinstance(choices2, list) and isinstance(choices2[0], dict):
                    finish_reason = choices2[0].get("finish_reason") or choices2[0].get("native_finish_reason")
                else:
                    finish_reason = None
            except Exception:
                finish_reason = None

        # final text is combined (if continuation happened) or initial text
        final_text = combined or text or ""
        if not final_text:
            raise RuntimeError("Empty response from OpenRouter (possible quota or provider issue).")
        # strip PART/CONTINUE markers if present
        part_match = _PART_RE.search(final_text)
        cont_match = _CONTINUE_TOKEN_RE.search(final_text)
        cut_index = None
        if part_match:
            cut_index = part_match.start()
        elif cont_match:
            cut_index = cont_match.start()
        if cut_index:
            final_text = final_text[:cut_index].strip()

        # coerce to 4 paragraphs best-effort
        paragraphs = [p.strip() for p in re.split(r"\n\s*\n", final_text) if p.strip()]
        if len(paragraphs) < 4:
            sentences = re.split(r'(?<=[.!?])\s+', final_text.strip())
            if len(sentences) >= 4:
                n = len(sentences)
                q = n // 4
                new_pars = []
                idx = 0
                for i in range(4):
                    take = q + (1 if i < (n % 4) else 0)
                    seg = " ".join(sentences[idx: idx + take]).strip()
                    if seg:
                        new_pars.append(seg)
                    idx += take
                paragraphs = new_pars
        exec_summary = "\n\n".join(paragraphs[:4]) if paragraphs else final_text.strip()
        return exec_summary.strip()
    except Exception as e:
        raise RuntimeError(f"Executive summary generation failed (OpenRouter): {e}") from e


# ----------------------
# PDF assembly (unchanged functionality, but coerce None -> "Not specified")
# and: Removed exon length & patient length columns from "Molecular Profiling" table per request.
# ----------------------
def _register_bell_mt(bell_mt_ttf: Optional[str] = None) -> str:
    if not REPORTLAB_AVAILABLE:
        return "Helvetica"
    if bell_mt_ttf and os.path.isfile(bell_mt_ttf):
        try:
            pdfmetrics.registerFont(TTFont("BellMT", bell_mt_ttf))
            return "BellMT"
        except Exception:
            pass
    return "Helvetica"


def generate_pdf_report(analysis: Dict, patient_info: Dict, doctor_notes: Optional[str],
                        exec_summary: str, out_path: str = "genomic_report.pdf", bell_mt_ttf: Optional[str] = None,
                        page_size=A4) -> str:
    if not REPORTLAB_AVAILABLE:
        raise RuntimeError("reportlab not available. Install: pip install reportlab")

    # Normalize values that may be None
    def _safe_str(v: Any, fallback="Not specified"):
        return str(v) if (v is not None and str(v).strip() != "") else fallback

    font_name = _register_bell_mt(bell_mt_ttf)
    doc = SimpleDocTemplate(out_path, pagesize=page_size,
                            leftMargin=20*mm, rightMargin=20*mm, topMargin=20*mm, bottomMargin=20*mm)
    styles = getSampleStyleSheet()
    normal = styles["Normal"]
    normal.fontName = font_name
    normal.fontSize = 10
    normal.leading = 13

    heading = ParagraphStyle('Heading', parent=styles['Heading1'])
    heading.fontName = font_name
    heading.fontSize = 18
    heading.leading = 22
    heading.alignment = 1  # center

    mono_style = ParagraphStyle('Mono', parent=styles.get('Code', styles['Normal']))
    mono_style.fontName = "Courier"
    mono_style.fontSize = 8
    mono_style.leading = 10

    small = ParagraphStyle('Small', parent=styles['Normal'])
    small.fontName = font_name
    small.fontSize = 9
    small.leading = 11

    story = []
    title_text = "<b>Patient Molecular Report</b>"
    story.append(Paragraph(title_text, heading))
    story.append(Spacer(1, 6))

    # Patient Info Table
    age = _safe_str(patient_info.get("age") if patient_info else None)
    sex = _safe_str(patient_info.get("sex") if patient_info else None)
    cancer = _safe_str(patient_info.get("cancer_type") if patient_info else None)
    notes = _safe_str(doctor_notes, fallback="Not specified")

    pat_table_data = [
        ["Field", "Value"],
        ["Age", age],
        ["Sex", sex],
        ["Cancer type", cancer],
        ["Doctor notes", notes],
    ]
    pat_table = Table(pat_table_data, colWidths=[60*mm, 80*mm])
    pat_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#f3f4f6")),
        ('TEXTCOLOR', (0,0), (-1,0), colors.black),
        ('FONTNAME', (0,0), (-1,-1), font_name),
        ('INNERGRID', (0,0), (-1,-1), 0.25, colors.grey),
        ('BOX', (0,0), (-1,-1), 0.5, colors.grey),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('ALIGN',(0,0),(-1,0),'CENTER')
    ]))
    story.append(Paragraph("<b>Patient Info</b>", small))
    story.append(Spacer(1, 4))
    story.append(pat_table)
    story.append(Spacer(1, 12))

    # Molecular Profiling table
    # NOTE: Exon length (ref) and Patient length columns removed as requested.
    gene = _safe_str(analysis.get("gene") if analysis else None, fallback="N/A")
    ref = _safe_str(analysis.get("ref_accession") if analysis else None, fallback="N/A")
    exon = _safe_str(analysis.get("exon") if analysis else None, fallback="N/A")
    variants_list = analysis.get("variants_hgvs") or []
    variants_text = ", ".join(variants_list) if variants_list else "None"
    mol_table_data = [
        ["Gene", "RefSeq mRNA", "Exon", "Variants"],
        [gene, ref, exon, variants_text]
    ]
    mol_table = Table(mol_table_data, colWidths=[30*mm, 70*mm, 20*mm, 55*mm])
    mol_table.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#f8fafc")),
        ('FONTNAME', (0,0), (-1,-1), font_name),
        ('INNERGRID', (0,0), (-1,-1), 0.25, colors.grey),
        ('BOX', (0,0), (-1,-1), 0.5, colors.grey),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('ALIGN',(0,0),(-1,0),'CENTER')
    ]))
    story.append(Paragraph("<b>Molecular Profiling</b>", small))
    story.append(Spacer(1, 4))
    story.append(mol_table)
    story.append(Spacer(1, 12))

    # Alignment excerpt (monospace)
    aligned_ref = (analysis.get("aligned_ref") or "").replace(" ", "") if analysis else ""
    aligned_obs = (analysis.get("aligned_obs") or "").replace(" ", "") if analysis else ""
    if aligned_ref or aligned_obs:
        story.append(Paragraph("<b>Alignment (reference vs patient)</b>", small))
        story.append(Spacer(1, 6))

        def chunk_text(s: str, width=100):
            if not s:
                return []
            return [s[i:i+width] for i in range(0, len(s), width)]

        ref_lines = chunk_text(aligned_ref, width=100)
        obs_lines = chunk_text(aligned_obs, width=100)
        max_lines = max(len(ref_lines), len(obs_lines))
        for i in range(max_lines):
            rline = ref_lines[i] if i < len(ref_lines) else ""
            oline = obs_lines[i] if i < len(obs_lines) else ""
            story.append(Paragraph(f"<font face='Courier'>{rline}</font>", mono_style))
            story.append(Paragraph(f"<font face='Courier'>{oline}</font>", mono_style))
        story.append(Spacer(1, 12))

    # Executive Summary (AI)
    story.append(Paragraph("<b>Executive Summary</b>", small))
    story.append(Spacer(1, 6))
    exec_pars = [p.strip() for p in re.split(r"\n\s*\n", exec_summary) if p.strip()] if exec_summary else []
    if not exec_pars:
        exec_pars = [exec_summary.strip()] if exec_summary else ["Executive summary not available."]
    for p in exec_pars:
        story.append(Paragraph(p, normal))
        story.append(Spacer(1, 6))

    story.append(Spacer(1, 12))
    footer = Paragraph("Report generated by Mutation Tool — AI Executive Summary (OpenRouter). Local tables generated from user input.", small)
    story.append(footer)

    doc.build(story)
    return os.path.abspath(out_path)


# ----------------------
# Combined wrapper
# ----------------------
def create_full_report_and_save(analysis: Dict, patient_info: Dict, doctor_notes: Optional[str] = None,
                                out_path: str = None, use_ai: bool = True,
                                bell_mt_ttf: Optional[str] = None, ai_kwargs: Optional[Dict] = None,
                                exec_summary: Optional[str] = None) -> str:
    """
    Create the full report PDF and save locally.

    exec_summary (optional): if provided, this text will be used as the executive summary
        (no AI call). If not provided and use_ai==True, the function will call generate_executive_summary.
    """
    if out_path is None:
        fd, out_path = tempfile.mkstemp(suffix=".pdf")
        os.close(fd)

    if ai_kwargs is None:
        ai_kwargs = {}

    # If an exec_summary string was provided by the caller (e.g., app preview), use it and skip AI call.
    exec_summary_text = None
    if exec_summary and isinstance(exec_summary, str) and exec_summary.strip():
        exec_summary_text = exec_summary.strip()
    elif use_ai:
        try:
            exec_summary_text = generate_executive_summary(analysis, patient_info, doctor_notes, **ai_kwargs)
        except Exception as e:
            exec_summary_text = f"Executive summary generation unavailable: {e}"
    else:
        exec_summary_text = "Executive summary generation disabled (use_ai=False)."

    # assemble PDF locally using the exec_summary_text
    pdf_path = generate_pdf_report(analysis, patient_info, doctor_notes, exec_summary_text, out_path=out_path, bell_mt_ttf=bell_mt_ttf)
    return pdf_path


if __name__ == "__main__":
    # quick local test (no AI by default)
    sample_analysis = {
        "gene": "KRAS",
        "ref_accession": "NM_001369786.1",
        "exon": 2,
        "exon_ref_length": 122,
        "patient_length": 122,
        "variants_hgvs": ["c.34G>T (p.G12C)"],
        "aligned_ref": "GCCTGCTGAAAATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAG",
        "aligned_obs": "GCCTGCTGAAAATGACTGAATATAAACTTGTGGTAGTTGGAGCTTGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAG"
    }
    sample_patient = {"age": "45", "sex": "MALE", "cancer_type": "lung adenocarcinoma"}
    sample_notes = "Smoker; chronic cough."
    try:
        out = create_full_report_and_save(sample_analysis, sample_patient, sample_notes, out_path="test_genomic_report.pdf", use_ai=False, exec_summary="Sample executive summary (test).")
        print("PDF written to:", out)
    except Exception as e:
        print("Failed to create PDF:", e)
