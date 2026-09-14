"""下载候选 GSE 的 suppl RAW.tar 到 raw_ext/。"""
import urllib.request, os, sys

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CAND = {"GSE199318":("GSE199nnn",None),"GSE95426":("GSE95nnn",None),
        "GSE165798":("GSE165nnn",None),"GSE303920":("GSE303nnn",None)}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "raw_ext")
os.makedirs(OUT, exist_ok=True)

def main(only=None):
    for gse,(folder,_) in CAND.items():
        if only and gse not in only:
            continue
        url = f"{BASE}/{folder}/{gse}/suppl/{gse}_RAW.tar"
        dst = os.path.join(OUT, f"{gse}_RAW.tar")
        print(f"GET {gse} <- {url}", file=sys.stderr)
        req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=600) as r:
            data = r.read()
        open(dst,"wb").write(data)
        print(f"  -> {len(data)} bytes", file=sys.stderr)

if __name__ == "__main__":
    only = sys.argv[1:] or None
    main(only)
