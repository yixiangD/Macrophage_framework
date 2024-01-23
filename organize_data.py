import subprocess
import os

for root, path, file in os.walk("./data"):
    root, path, file = root, path, file
    break

data10x = []
for p in path:
    for rr, pp, ff in os.walk(os.path.join(root, p)):
        sc_files = ff
        sc_folders = pp
        break
    if len(pp) > 0:
        print(pp)
    for f in ff:
        if "barcode" not in f and "matrix" not in f \
                and "gene" not in f and "feature" not in f:
            if "dge" not in f:
                if "RData.gz" in f:
                    # unzip RData.gz
                    tarpath = os.path.join(rr, f)
                    subprocess.call(f"gzip -dv {tarpath}", shell=True)
                else:
                    # read count tables
                    pass
            else:
                # read dge files
                temp = f.replace(".dge.txt.gz", "")
                # print(temp.split("_", 1)[1])
        else:
            # standard 10x data
            # default: barcodes, genes, matrix
            # need to replace "features" to "genes"
            data10x.append([p] + f.rsplit("_", 1))
import pandas as pd
df10x = pd.DataFrame(data10x, columns=["folder", "group", "file"])
df10x[["sample", "type"]] = df10x["group"].str.split("_", n = 1, expand=True)
print(df10x)

df10x_count = df10x.groupby("folder")["group"].agg(list).reset_index()
df10x_count["count"] = df10x_count ["group"].apply(len)
#print(df10x_count)

df10x_count = df10x.groupby(["folder", "sample"])["type"].agg(list).reset_index()
df10x_count["count"] = df10x_count ["type"].apply(len)
print(df10x_count)
