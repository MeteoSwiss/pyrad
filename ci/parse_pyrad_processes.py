# %%
import pandas as pd
from pathlib import Path
import os
import inspect
import yaml
import re
from pyrad import proc

mainpath = Path(__file__).resolve().parent.parent
procpath = Path(mainpath, 'src', 'pyrad_proc', 'pyrad', 'proc')
OUT_DIRECTORY = str(Path(mainpath, 'doc', 'source', 'overview'))
DOC_PATH = 'https://meteoswiss.github.io/pyrad/_modules/pyrad'

def funcpath_to_docpath(funcpath):
    funcpath = funcpath.split('pyrad')[-1]
    return DOC_PATH + funcpath.replace('.py', '.html')

def dict_to_restructured_text(yaml_data):
    rst_output = []
    rst_output.append('List of pyrad datasets')
    rst_output.append('==============================\n')

    for key, value in yaml_data.items():
        rst_output.append(f"{key}")
        rst_output.append("-----------------------------")

        for key2, value2 in yaml_data[key].items():
            if 'description' not in value2.keys():
                continue
            rst_output.append(f"{key2}")
            rst_output.append('""""""""""""""""""""""""""""""')
            rst_output.append('description')
            rst_output.append('   ' + value2['description'] + f'\n `[Source] <{value2["link"]}>`_' )
            rst_output.append('parameters')
            params = value2['parameters']
            for key3, value3 in params.items():
                rst_output.append('   ' + key3)
                rst_output.append('       ' + value3)
            rst_output.append('')
            # rst_output.append(f"\n\n{value2['parameters']}\n\n")
    return '\n'.join(rst_output)


def parse_string_to_dict(input_string):
    lines = input_string.strip().split('\n')
    result_dict = {}
    current_key = None

    for line in lines:
        if line.startswith(' '):
            line = line[4:]
        if line.startswith(" "):
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
    description = ''
    read_desc = True
    for line in lines:
        if 'Parameters' in line:
            read_desc = False
        if read_desc:
            description += line + ' '
        if ('Accepted Configuration Keywords:'.lower() in line.lower()):
            start_reading_attr = True
            continue
        if 'radar_list' in line:
            start_reading_attr = False

        if start_reading_attr:
            attributes += line + '\n'

    dic = {'description': " ".join(description.split()),
           'parameters': parse_string_to_dict(attributes)}
    return dic


# Match function names to product types
all_processes = {}
mapping_func_to_prodname = {}
file_path = Path(procpath, 'process_aux.py')
with open(file_path, 'r') as f:
    content = f.readlines()
    for line in content:
        if 'format output' in line:
            processtype = line.split('format output')[
                0].strip().replace("'", "")
            all_processes[processtype] = {}
        match = re.findall("^'[A-Z0-9_]*'\\s*:\\d*", line.strip())
        if len(match):
            process = match[0].split(':')[0].replace("'", "")
            function = line.split(':')[1].strip()
            all_processes[processtype][process] = function

for processtype in all_processes:
    for process in all_processes[processtype]:
        fun = getattr(proc, all_processes[processtype][process])
        funcname = all_processes[processtype][process]
        docstr = inspect.getdoc(fun)
        parameters = process_docstring(docstr)
        parameters['link'] = funcpath_to_docpath(inspect.getfile(fun)) +'#' + funcname 
        all_processes[processtype][process] = parameters

# Convert dict data to reStructuredText format
rst_content = dict_to_restructured_text(all_processes)
fname = Path(OUT_DIRECTORY, 'list_process.rst')
with open(fname, 'w') as f:
    f.write(rst_content)


# %%
