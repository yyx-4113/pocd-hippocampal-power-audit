"""
核验候选 bulk 数据集 (E-utilities, HTTPS)。
目标: 判定 GSE199318 / GSE95426 / GSE165798 / GSE303920 / GSE316433 / GSE234493
是否可用作"海马 bulk POCD 元分析"的扩展。
"""
import urllib.request, urllib.parse, json, re, time, sys

EUTIL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACCS = ["GSE199318", "GSE95426", "GSE165798", "GSE303920", "GSE316433", "GSE234493"]

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

def esearch(acc):
    url = f"{EUTIL}/esearch.fcgi?db=gds&term={acc}[accn]&retmode=json"
    d = json.loads(fetch(url))
    ids = d.get("esearchresult", {}).get("idlist", [])
    return ids[0] if ids else None

def esummary(uid):
    url = f"{EUTIL}/esummary.fcgi?db=gds&id={uid}&retmode=json"
    d = json.loads(fetch(url))
    return d.get("result", {}).get(uid, {})

def main():
    out = []
    for acc in ACCS:
        print(f"== {acc} ==", file=sys.stderr)
        try:
            uid = esearch(acc)
            if not uid:
                out.append({"acc": acc, "uid": None, "error": "no uid"})
                continue
            s = esummary(uid)
            rec = {
                "acc": acc,
                "uid": uid,
                "title": s.get("title", ""),
                "summary": (s.get("summary", "") or "")[:600],
                "organism": s.get("organism", ""),
                "gpl": s.get("GPL", ""),
                "platform": s.get("platform", ""),
                "samples": s.get("samples", ""),
                "sample_count": s.get("sample_count", ""),
                "gse": s.get("GSE", acc),
                "submission": s.get("submission", ""),
                "pubdate": s.get("PDAT", ""),
            }
            out.append(rec)
        except Exception as e:
            out.append({"acc": acc, "error": repr(e)})
        time.sleep(0.4)
    print(json.dumps(out, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main()
