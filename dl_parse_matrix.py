"""下载 4 个候选 GSE series matrix 并解析分组/平台。"""
import urllib.request, gzip, os, io, re

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CAND = {"GSE199318":"GSE199nnn","GSE95426":"GSE95nnn","GSE165798":"GSE165nnn","GSE303920":"GSE303nnn"}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix_tmp")
os.makedirs(OUT, exist_ok=True)

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()

def parse(text):
    lines = text.splitlines()
    meta = {}
    cols = []
    data_start = None
    for i, ln in enumerate(lines):
        if ln.startswith("!") or ln.startswith("^"):
            if ln.startswith("!Sample_title"):
                meta["title"] = ln.split("\t")[1:]
            elif ln.startswith("!Sample_characteristics_ch1"):
                meta["char"] = ln.split("\t")[1:]
            elif ln.startswith("!Sample_source_name_ch1"):
                meta["source"] = ln.split("\t")[1:]
            elif ln.startswith("!Sample_geo_accession"):
                meta["gsm"] = ln.split("\t")[1:]
            elif ln.startswith("!Platform_title") or ln.startswith("!Platform_organism"):
                meta.setdefault("platform", []).append(ln[1:])
        else:
            if data_start is None and i>0 and lines[i-1].startswith("!series_matrix_table_begin"):
                data_start = i
            if ln.startswith("!series_matrix_table_end"):
                data_start = None
    return meta

def main():
    for gse, folder in CAND.items():
        url = f"{BASE}/{folder}/{gse}/matrix/{gse}_series_matrix.txt.gz"
        raw = get(url)
        open(os.path.join(OUT, f"{gse}.series_matrix.txt.gz"), "wb").write(raw)
        text = gzip.decompress(raw).decode("utf-8", "replace")
        m = parse(text)
        print(f"\n########## {gse} ##########")
        print("platform:", m.get("platform"))
        print("n_samples:", len(m.get("title", [])))
        for j in range(len(m.get("title", []))):
            t = m.get("title", [""]*99)[j]
            c = m.get("char", [""]*99)[j] if "char" in m else ""
            s = m.get("source", [""]*99)[j] if "source" in m else ""
            print(f"  [{j}] title={t!r} | char={c!r} | src={s!r}")

if __name__ == "__main__":
    main()
