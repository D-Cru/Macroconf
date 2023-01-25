# Convenience function to rename hashes, e.g. if new parameters are added.
import os

old = [
    "old_hash_1",
    "old_hash_2",
]


new = [
    "new_hash_1",
    "new_hash_2",
]


for o, n in zip(old, new):
    print(o, n)
    # snakemake-config.yaml reports/paper_figures.csv
    output = f"sed -i 's/{o}/{n}/g' reports/paper_figures.csv"
    # output=f"find . -depth -exec rename -v 's/{o}/{n}/g' {{}} + &"
    # print(output)
    os.system(output)
    # os.system("echo")
