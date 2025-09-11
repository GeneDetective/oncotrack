________________________________________
# OncoTrack — Detect. Decode. Decide.

OncoTrack is a web app that fetches a canonical RefSeq transcript from the GRCh38 dataset, extracts an exon region, and compares a patient DNA exon sequence (FASTA or raw A/C/G/T/N) to detect variants (mutation). The app produces HGVS-like output and generates physician-friendly summaries using an LLM back-end (OpenRouter / GPT-OSS-20B).

> **Demo:** → https://oncotrack.onrender.com/
> **Note:** This repository uses synthetic / literature-based test sequences for demonstration. No protected health information or real patient data is included.

## Quick highlights

- Detects SNVs and simple indels inside a chosen exon.
- Uses **NCBI RefSeq** transcripts for canonical exon definitions.
- Uses **OpenRouter (gpt-oss-20b)** for report generation (requires an API key).
- Example validated hotspots included: KRAS G12C/G12D/G13D, NRAS Q61R, BRAF V600E, PIK3CA H1047R (see `tests/`).

---

## Requirements
- Python 3.9+ (3.10/3.11 recommended)
- Git
- Internet access (for NCBI queries and OpenRouter inference)
- Flask 3.0.3
- Biopython 1.85
- Requests 2.32.5
- Reportlab 4.2.5


---

## Setup
### 1. Clone the repo
```
git clone https://github.com/GeneDetective/oncotrack.git
```
```
cd oncotrack
code .
```
2. Create a virtual environment & activate it
Linux / macOS:
```
python -m venv venv
source venv/bin/activate
```
Windows (CMD):
```
python -m venv venv
venv\Scripts\activate
```
Windows (PowerShell):
```
python -m venv venv
.\venv\Scripts\Activate.ps1
```
3. Install dependencies
```
pip install --upgrade pip
pip install -r requirements.txt
```
4. Set environment variables

Environment variables required by the app:
```
•	NCBI_ENTREZ_EMAIL — your email (required by NCBI Entrez).
•	NCBI_ENTREZ_API_KEY — optional but recommended (NCBI API key).
•	OPENROUTER_API_KEY — OpenRouter API key for LLM inference.
•	OPENROUTER_MODEL — which model to use: openai/gpt-oss-20b:free
```
Linux / macOS (bash/zsh):
```
export NCBI_ENTREZ_EMAIL="your.email@example.com"
export NCBI_ENTREZ_API_KEY="your_ncbi_key_or_blank"
export OPENROUTER_API_KEY="your_openrouter_key"
export OPENROUTER_MODEL="openai/gpt-oss-20b:free"
```
Windows (CMD) using setx (persists across sessions):
```
setx NCBI_ENTREZ_EMAIL "your.email@example.com"
setx NCBI_ENTREZ_API_KEY "your_ncbi_key"
setx OPENROUTER_API_KEY "your_openrouter_key"
setx OPENROUTER_MODEL "openai/gpt-oss-20b:free"
```
After setx, you need to open a new terminal window for variables to be available.

Windows PowerShell (persistent for current user):
```
setx NCBI_ENTREZ_EMAIL "your.email@example.com"
setx NCBI_ENTREZ_API_KEY "your_ncbi_key"
setx OPENROUTER_API_KEY "your_openrouter_key"
setx OPENROUTER_MODEL "openai/gpt-oss-20b:free"
```
________________________________________
How to get the API keys 
**OpenRouter (for LLM)**
1.	Go to https://openrouter.ai/ and create an account.
2.	Visit their Models/APIs page (select openai/gpt-oss-20b:free).
3.	Create an API key (name it e.g. oncotrack) and copy it.
4.	Set OPENROUTER_API_KEY in environment variables as shown above.
5.	Set OPENROUTER_MODEL to openai/gpt-oss-20b:free.

**NCBI Entrez**
1.	Create or sign into an NCBI account: https://account.ncbi.nlm.nih.gov/
2.	Go to Settings → API key management and create a key.
3.	Set NCBI_ENTREZ_EMAIL (your email) and NCBI_ENTREZ_API_KEY environment variables.
________________________________________
Run locally
# make sure venv is activated and env vars are set
```python app.py```
Open your browser at http://127.0.0.1:5000 (or the port printed by the server).
```
Example input (paste into the form):
Gene: KRAS
Exon: 2
>patient_KRAS_exon2_G12C
GCCTGCTGAAAATGACTGAATATAAACTTGTGGTAGTTGGAGCTTGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGACGAATATGATCCAACAATAGAG
Age: 64
Sex: Male
Cancer Type: Lung adenocarcinoma
Doctor's notes: Smoker, persistent cough for 3 months.
```
**CLICK ON ANALYZE**
```
Expected output: c.34G>T (p.G12C)
```
Then;
**CLICK ON REPORT**
```
Executive Detailed Report Generation  
```
*You can try it out with other mock patient data from: https://drive.google.com/file/d/1uE37nf-VNRq3HOw3DbLpmhW2PYTxliyR/view?usp=sharing
________________________________________
# oncotrack
