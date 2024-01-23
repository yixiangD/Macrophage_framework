import subprocess
import os

for root, path, file in os.walk("./data"):
    root, path, file = root, path, file
    break

for p in path:
    for rr, pp, ff in os.walk(os.path.join(root, p)):
        sc_files = ff
        sc_folders = pp
        break
    if len(pp) > 0:
        print(pp)
    for f in ff:
        if "barcode" not in f and "matrix" not in f and "gene" not in f and "feature" not in f:
            if "dge" not in f:
                print(rr)
                print(f)
                if "RData.gz" in f:
                    # unzip RData.gz
                    tarpath = os.path.join(rr, f)
                    subprocess.call(f"gzip -dv {tarpath}", shell=True)
                else:
                    # read count tables
                    pass
