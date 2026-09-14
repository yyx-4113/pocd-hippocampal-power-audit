"""下载候选 GSE 的 series matrix 到本地, 用于确认分组与平台。"""
import urllib.request, os, sys

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CAND = {
    "GSE199318": "GSE199nnn",
    "GSE95426":  "GSE95nnn",
    "GSE165798": "GSE165nnn",
    "GSE303920": "GSE303nnn",
}
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix_tmp")
os.makedirs(OUTDIR, exist_ok=True)

def main():
    for gse, folder in CAND.items():
        url = f"{BASE}/{folder}/{gse}/matrix/{gse}.series_matrix.txt.gz"
        dst = os.path.join(OUTDIR, f"{gse}.series_matrix.txt.gz")
        print(f"GET {url}", file=sys.stderr)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            open(dst, "wb").write(data)
            print(f"  saved {len(data)} bytes -> {dst}")
        except Exception as e:
            print("  ERROR:", repr(e))

if __name__ == "__main__":
    main()
