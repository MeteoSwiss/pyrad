#%%
import pandas as pd
from pathlib import Path
import os
import inspect
import yaml
import re
from pyrad import proc

mainpath = Path(__file__).resolve().parent.parent
procpath = Path(mainpath, 'src', 'pyrad_proc', 'pyrad', 'proc')
OUT_DIRECTORY = str(Path(mainpath, 'doc', 'source', 'overview', 'mappings'))
            


def parse_string_to_dict(input_string):
    lines = input_string.strip().split('\n')
    result_dict = {}
    current_key = None

    for line in lines:
        if line.startswith(' '):
            line = line[4:]
        if line.startswith(" ") and ':' not in line:
            if current_key is not None:
                result_dict[current_key] += line
        else:
            current_key = line
            if current_key not in result_dict:
                result_dict[current_key] = ''
    
    for k in result_dict:
        result_dict[k] = result_dict[k].strip().strip('\n')
        result_dict[k] = " ".join(result_dict[k].split())

    return result_dict
    
def process_docstring(docstr):
    start_reading_attr = False
    lines = docstr.split('\n')

    attributes = ''
    for line in lines:
        if ('Accepted Configuration Keywords:'.lower() in line.lower()):
            start_reading_attr = True
            continue
        if 'radar_list' in line:
            start_reading_attr = False

        if start_reading_attr:
            attributes += line +'\n'
    return parse_string_to_dict(attributes)

# Match function names to product types
all_processes = {}
mapping_func_to_prodname = {}
file_path = Path(procpath, 'process_aux.py')
with open(file_path, 'r') as f:
    content = f.readlines()
    for line in content:
        if 'format output' in line:
            processtype = line.split('format output')[0].strip().replace("'","")
            all_processes[processtype] = {}
        match = re.findall("^'[A-Z0-9_]*'\s*:\d*", line.strip())
        if len(match):
            process = match[0].split(':')[0].replace("'","")
            function = line.split(':')[1].strip()
            all_processes[processtype][process] = function


for processtype in all_processes:
    for process in all_processes[processtype]:
        fun = getattr(proc, all_processes[processtype][process])
        docstr = inspect.getdoc(fun)
        para.0meters = process_docstring(docstr)


    
    
# %%
