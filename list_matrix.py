import urllib.request, re, sys
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CAND = {"GSE199318":"GSE199nnn","GSE95426":"GSE95nnn","GSE165798":"GSE165nnn","GSE303920":"GSE303nnn"}
def list_dir(url):
    req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8","replace")
    return [l for l in re.findall(r'href="([^"]+)"', html) if l not in ("","../")]
for gse, folder in CAND.items():
    for sub in ("matrix/", "miniml/", "soft/"):
        url = f"{BASE}/{folder}/{gse}/{sub}"
        try:
            ls = list_dir(url)
            print(f"{gse}/{sub}: {ls}")
        except Exception as e:
            print(f"{gse}/{sub}: ERR {e}")
