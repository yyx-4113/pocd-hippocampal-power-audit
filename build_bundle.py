# -*- coding: utf-8 -*-
"""把 geo_meta 工作目录组装成可公开的复现包（GitHub 仓库内容）。

设计原则
--------
1. **不改动任何已验证脚本**：脚本内部用 os.path.dirname(__file__) 解析路径，
   所以复现包必须**平铺**（仓库根 = 原工作目录）。这样 `python step5_sc_analysis.py`
   在仓库根直接跑通，输出落在 ./sc_out/，与归档结果一致。
2. **排除可再下载的大文件与大中间体**（GEO 原始数据、h5ad），README 里写明取回方式。
3. **排除稿件侧文档**（投稿信/清单/审稿报告）：仓库只放代码 + 结果，稿件另行投稿。

产物：<工作区根>/pocd-hippocampal-power-audit/
用法：python build_bundle.py [--dest 目标目录]
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEST = os.path.join(os.path.dirname(HERE), "pocd-hippocampal-power-audit")
if "--dest" in sys.argv:
    DEST = sys.argv[sys.argv.index("--dest") + 1]

# 目录级排除（体积大 / 可再下载 / 临时）
SKIP_DIRS = {
    "__pycache__", "raw_ext", "GSE267933_RAW",
    ".git", "sc_out_all", "pocd-hippocampal-power-audit",
}
# 文件级排除（精确名）
SKIP_FILES = {
    "_before_step5d.csv",
    "_audit.txt",                 # 命令行重定向产物；重跑 verify_manuscript_numbers.py 即得
    "GSE267933.h5ad",             # 549 MB，GEO 可下载，且由 step1b 重建
    "GSE267933_RAW.tar",          # 134 MB，GEO 原始 tar
    # 稿件侧工具（生成/校验投稿 docx 与 .md），不属于分析复现链
    "build_docx.py",
    "verify_docx_v3.py",
    "verify_md_v3.py",
}
# sc_out 内排除（715 MB 中间体）
SKIP_SCO = {"GSE267933_processed.h5ad"}
# 稿件侧文档（不随仓库公开）
SKIP_PREFIX = ("路线A+C", "审稿报告", "manuscript_AA.html", "Step5c_结论校准")
SKIP_SUFFIX = (".err",)          # 空错误日志


def keep(name):
    if name in SKIP_FILES:
        return False
    if name.startswith(SKIP_PREFIX):
        return False
    if name.endswith(SKIP_SUFFIX):
        return False
    return True


def main():
    if os.path.isdir(os.path.join(DEST, ".git")):
        print("[warn] 目标已有 .git，仅覆盖文件，不重建仓库")
    os.makedirs(DEST, exist_ok=True)
    n, skipped = 0, []

    for root, dirs, files in os.walk(HERE):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        rel = os.path.relpath(root, HERE)
        if rel == ".":
            rel = ""
        # 不在包内嵌套包
        if DEST.startswith(root + os.sep):
            continue
        for d in list(dirs):
            if os.path.join(root, d) == DEST:
                dirs.remove(d)
        out_root = os.path.join(DEST, rel) if rel else DEST
        for fn in files:
            if not keep(fn):
                skipped.append(os.path.join(rel, fn))
                continue
            if rel == "sc_out" and fn in SKIP_SCO:
                skipped.append(os.path.join(rel, fn))
                continue
            os.makedirs(out_root, exist_ok=True)
            shutil.copy2(os.path.join(root, fn), os.path.join(out_root, fn))
            n += 1

    print(f"[done] 复制 {n} 个文件 -> {DEST}")
    print(f"[skip] {len(skipped)} 项：")
    for s in sorted(skipped):
        print("   -", s)
    tot = 0
    for dp, dns, fs in os.walk(DEST):
        dns[:] = [d for d in dns if d != ".git"]
        tot += sum(os.path.getsize(os.path.join(dp, f)) for f in fs)
    print(f"[size] 复现包体积 {tot / 1e6:.1f} MB（不含 .git）")


if __name__ == "__main__":
    main()
