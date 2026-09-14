"""列出 GEO FTP 系列目录, 确认可下载的处理矩阵文件。"""
import urllib.request, re, sys

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
CAND = {
    "GSE199318": "GSE199nnn",
    "GSE95426":  "GSE95nnn",
    "GSE165798": "GSE165nnn",
    "GSE303920": "GSE303nnn",
}

def list_dir(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    links = re.findall(r'href="([^"]+)"', html)
    return [l for l in links if not l.endswith("../") and l not in ("",)]

def main():
    for gse, folder in CAND.items():
        url = f"{BASE}/{folder}/{gse}/"
        print(f"\n===== {gse}  ({url}) =====")
        try:
            ls = list_dir(url)
            for l in ls:
                print(" ", l)
            # 也列 suppl
            sup = f"{url}suppl/"
            try:
                ls2 = list_dir(sup)
                print(f"  -- suppl/ ({len(ls2)} files) --")
                for l in ls2[:40]:
                    print("    ", l)
            except Exception as e:
                print("  (no suppl or err:", repr(e), ")")
        except Exception as e:
            print("  ERROR:", repr(e))

if __name__ == "__main__":
    main()
