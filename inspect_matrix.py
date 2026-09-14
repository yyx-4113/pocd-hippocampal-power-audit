"""检查各 series matrix 的行标识与量纲。"""
import gzip, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix_tmp")
for gse in ["GSE199318","GSE95426","GSE165798","GSE303920"]:
    p = os.path.join(OUT, f"{gse}.series_matrix.txt.gz")
    text = gzip.decompress(open(p,"rb").read()).decode("utf-8","replace")
    lines = text.splitlines()
    # 找表头行 (!series_matrix_table_begin 的下一行是列名; 再下一行是数据)
    begin = [i for i,l in enumerate(lines) if l.startswith("!series_matrix_table_begin")]
    hdr_i = begin[0]+1
    header = lines[hdr_i].split("\t")
    print(f"\n########## {gse} ##########")
    print("n_cols:", len(header))
    print("col0:", repr(header[0]), "| samples:", header[1:4], "...", header[-2:])
    # 数据行
    for k in range(hdr_i+1, hdr_i+6):
        row = lines[k].split("\t")
        vals = row[1:]
        # 简单统计
        try:
            fv = [float(v) for v in vals]
            desc = f"min={min(fv):.2f} max={max(fv):.2f} mean={sum(fv)/len(fv):.2f}"
        except:
            desc = "non-numeric"
        print(f"  rowid={row[0][:40]!r}  {desc}")
    # 平台 id 行
    for l in lines:
        if l.startswith("!Platform_ID") or l.startswith("!platform_id"):
            print("platform_id_field:", l[:80]); break
