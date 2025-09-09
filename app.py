# app.py
from flask import Flask, render_template, request, redirect, url_for, flash
from sequence_utils import analyze_patient_vs_reference, AnalysisError
import os
import traceback

# Robust import logic:
# Prefer a generate_executive_summary if available; no PDF generation is used.
generate_executive_summary = None
generate_report_legacy = None
try:
    from ai_report import generate_executive_summary as _gen_exec  # type: ignore
    generate_executive_summary = _gen_exec
except Exception:
    try:
        # legacy fallback if available on your server
        from ai_report import generate_report as _gen_legacy  # type: ignore
        generate_report_legacy = _gen_legacy
    except Exception:
        pass

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET", "replace-with-a-random-secret")


@app.route("/", methods=["GET", "POST"])
def index():
    """
    Landing page: runs the sequence analysis and shows a preview/launchpad to create the AI report.
    Keeps the pipeline unchanged.
    """
    result = None
    error = None
    if request.method == "POST" and request.form.get("action") == "analyze":
        gene = request.form.get("gene", "").strip()
        exon = request.form.get("exon", "").strip()
        patient_seq = request.form.get("sequence", "").strip()
        age = request.form.get("age", "").strip()
        sex = request.form.get("sex", "").strip()
        cancer_type = request.form.get("cancer_type", "").strip()
        # capture doctor's notes so they persist into the result and can be passed on to report
        doctor_notes = request.form.get("doctor_notes", "").strip()

        try:
            analysis = analyze_patient_vs_reference(gene, exon, patient_seq)
            result = {
                "gene": analysis["gene"],
                "ref_accession": analysis.get("ref_accession", ""),
                "exon": analysis.get("exon", exon),
                "exon_ref_length": analysis.get("exon_ref_length"),
                "patient_length": analysis.get("patient_length"),
                "variants_hgvs": analysis.get("variants_hgvs", []),
                "aligned_ref": analysis.get("aligned_ref", ""),
                "aligned_obs": analysis.get("aligned_obs", ""),
                "patient_info": {"age": age, "sex": sex, "cancer_type": cancer_type},
                "note": analysis.get("note"),
                "doctor_notes": doctor_notes,
            }
        except AnalysisError as e:
            error = str(e)
        except Exception as e:
            error = f"Unexpected error: {e}\n{traceback.format_exc()}"

    return render_template("index.html", result=result, error=error)


@app.route("/generate_report", methods=["POST"])
def generate_report_route():
    """
    Build the in-browser AI report preview (no PDFs). We keep the preview-first behavior:
    - Collect analysis + patient info + notes from the form
    - Generate the AI executive summary (or fallback text)
    - Render templates/report_preview.html with tables + summary
    """
    try:
        # helper to read doctor notes (multiple possible fields)
        def get_field(*names):
            for n in names:
                v = request.form.get(n)
                if v is not None and str(v).strip() != "":
                    return v
            return ""

        gene = request.form.get("gene", "").strip()
        exon = request.form.get("exon", "").strip()
        ref_accession = request.form.get("ref_accession", "").strip()
        variants_hgvs_raw = request.form.get("variants_hgvs", "").strip()
        aligned_ref = request.form.get("aligned_ref", "").strip()
        aligned_obs = request.form.get("aligned_obs", "").strip()
        age = request.form.get("age", "").strip()
        sex = request.form.get("sex", "").strip()
        cancer_type = request.form.get("cancer_type", "").strip()
        # Prefer hidden_doctor_notes (populated on analyze), fallback to visible fields
        doctor_notes = get_field("hidden_doctor_notes", "doctor_notes", "notes", "doctorNotes")

        # Normalize variants input into a list (mirrors server logic you had)
        variants_hgvs = []
        if variants_hgvs_raw:
            for sep in ("\n", ";", ","):
                if sep in variants_hgvs_raw:
                    variants_hgvs = [v.strip() for v in variants_hgvs_raw.split(sep) if v.strip()]
                    break
            else:
                variants_hgvs = [variants_hgvs_raw.strip()]

        # Read exon lengths added as hidden fields in the form (from analyze step).
        # (We do not display these in the Molecular Profiling table per your earlier change.)
        exon_ref_length = request.form.get("exon_ref_length") or None
        patient_length = request.form.get("patient_length") or None

        analysis = {
            "gene": gene,
            "exon": exon,
            "ref_accession": ref_accession,
            "variants_hgvs": variants_hgvs,
            "aligned_ref": aligned_ref,
            "aligned_obs": aligned_obs,
            "exon_ref_length": exon_ref_length,
            "patient_length": patient_length,
        }
        patient_info = {"age": age or None, "sex": sex or None, "cancer_type": cancer_type or None}

        # ---------- Preview text (Executive Summary) ----------
        report_text = ""
        if generate_executive_summary is not None:
            try:
                # Generate the AI text (this is the preview content).
                # Use a slightly conservative token budget so longer inputs don't overflow.
                report_text = generate_executive_summary(
                    analysis, patient_info, doctor_notes or None, desired_output_tokens=650
                )
            except Exception as e:
                app.logger.exception("generate_executive_summary failed")
                # Preserve the detailed error message for transparency in the preview.
                report_text = f"Executive summary generation failed: {e}"
        elif generate_report_legacy is not None:
            try:
                report_text = generate_report_legacy(analysis, patient_info, doctor_notes or None)
            except Exception as e:
                app.logger.exception("generate_report (legacy) failed")
                report_text = f"Report generation failed: {e}"
        else:
            report_text = "AI report generation not available on this server."

        # ---------- Render the HTML-only preview (no PDF path at all) ----------
        return render_template(
            "report_preview.html",
            report_text=report_text,
            analysis=analysis,
            patient_info=patient_info,
            doctor_notes=doctor_notes or ""
        )

    except Exception as e:
        app.logger.exception("Failed to generate report.")
        flash(f"Failed to generate report: {e}")
        return redirect(url_for("index"))


# Retained utility; harmless if unused elsewhere.
def wrap_text(text: str, max_chars: int):
    words = text.split()
    if not words:
        return []
    lines = []
    cur = []
    cur_len = 0
    for w in words:
        additional = len(w) + (1 if cur else 0)
        if cur_len + additional > max_chars:
            lines.append(" ".join(cur))
            cur = [w]
            cur_len = len(w)
        else:
            cur.append(w)
            cur_len += additional
    if cur:
        lines.append(" ".join(cur))
    return lines


if __name__ == "__main__":
    app.run(debug=True)
