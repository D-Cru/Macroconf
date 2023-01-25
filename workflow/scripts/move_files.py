import os
import subprocess

# read file old.txt
with open("untracked.txt", "r") as f:
    old = f.readlines()

with open("untracked_destination.txt", "r") as f:
    new = f.readlines()


for o, n in zip(old, new):
    o = o.rstrip()
    n = n.rstrip()
    # print(o, n)
    # print(os.path.dirname(n))
    # get the file path from o
    output = f"mkdir -p {os.path.dirname(n)}"
    # print(output)
    # os.system(output)
    # output = f"mv -i {o} {n}"
    output = ["mv", "-i", o, n]
    # output = ['pwd']
    print(output)
    subprocess.call(output)
    # os.system(output)
#     # output=f"find . -depth -exec rename -v 's/{o}/{n}/g' {{}} + &"
#     # print(output)
#     os.system(output)
#     # os.system("echo")
