import subprocess
import os

for root, path, file in os.walk("./data"):
    root, path, file = root, path, file
    break
print(root, path, file)

for f in file:
    if f.endswith("RAW.tar"):
        fname = f.split("_")[0]
        fdir = os.path.join(root, fname)
        tarpath = os.path.join(root, f)

        if not os.path.exists(fdir):
            os.makedirs(fdir)
        subprocess.call(f"tar -C {fdir} -xzvf {tarpath}", shell=True)
