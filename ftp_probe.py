import ftplib, sys

def ls(remote_dir):
    try:
        ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov", timeout=30)
        ftp.login()
        print(f"### DIR: {remote_dir}")
        try:
            lines = ftp.nlst(remote_dir)
            for l in lines:
                print("  ", l)
        except Exception as e:
            print("  nlst error:", e)
        ftp.quit()
    except Exception as e:
        print(f"CONNECT ERROR for {remote_dir}: {e}")

for d in [
    "/geo/series/GSE215nnn/GSE215410/",
    "/geo/series/GSE215nnn/GSE215410/suppl/",
    "/geo/series/GSE215nnn/GSE215410/matrix/",
    "/geo/series/GSE178nnn/GSE178995/",
    "/geo/series/GSE178nnn/GSE178995/suppl/",
    "/geo/series/GSE178nnn/GSE178995/matrix/",
]:
    ls(d)
