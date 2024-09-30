#%%

from pathlib import Path
import inspect
import re
from pyrad import proc
import os

all_with_voltype = []
def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.readlines()
    cprod = None
    for line in content:
        if 'prdcfg["type"] ==' in line:
            cprod = line.split('==')[-1].strip()
        if 'voltype' in line:
            if cprod not in all_with_voltype:
                all_with_voltype.append(cprod)
            
mainpath = Path(__file__).resolve().parent.parent
prodpath = Path(mainpath, 'src', 'pyrad_proc', 'pyrad', 'prod')


for root, _, files in os.walk(prodpath):
    for file in files:
        if file.endswith(".py") and file.startswith("process"):
            file_path = os.path.join(root, file)
            print(file_path)
            process_file(file_path)
            
with open('all_with_voltype.txt','w') as f:
    for prod in all_with_voltype:
        print(prod)
        f.write(prod + '\n')
            