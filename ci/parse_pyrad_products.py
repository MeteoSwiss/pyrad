# %%
import pandas as pd
from pathlib import Path
import os 
import re

mainpath = Path(__file__).resolve().parent.parent
prodpath = Path(mainpath, 'src', 'pyrad_proc', 'pyrad', 'prod')
OUT_DIRECTORY = str(Path(mainpath, 'doc', 'source', 'overview'))
DOC_PATH = 'https://github.com/MeteoSwiss/pyrad/blob/master/src/pyrad_proc/pyrad'

def funcpath_to_docpath(funcpath):
    funcpath = funcpath.split('pyrad')[-1]
    return DOC_PATH + funcpath

def parameters_to_dict(params):
    dic = {}
    keys = re.findall('([a-zA-Z_]*\\s*:\\s*[a-zA-Z]*)', params)
    for i in range(len(keys)):
        idx_now = params.index(keys[i]) + len(keys[i])
        if i == len(keys) - 1:
            idx_next = len(params)
        else:
            idx_next = params.index(keys[i + 1])

        dic[keys[i].replace(':', '').strip()] = params[idx_now:idx_next]
    return dic


def dict_to_restructured_text(yaml_data):
    rst_output = []
    rst_output.append('List of pyrad products')
    rst_output.append('==============================\n')

    for datasettype, value in yaml_data.items():
        rst_output.append(f"{datasettype}")
        rst_output.append("-----------------------------")

        if 'all_products' in yaml_data[datasettype]:
            rst_output.append(".. note::")
            rst_output.append("   Supports all products of type " + \
                f"{yaml_data[datasettype]['all_products']}")
            rst_output.append('')

        for product, prodinfo in yaml_data[datasettype].items():
            if type(prodinfo) == str:
                continue
            if 'description' not in prodinfo:
                continue
            rst_output.append(f"{product}")
            rst_output.append('""""""""""""""""""""""""""""""')
            rst_output.append('description')
            rst_output.append('   ' + prodinfo['description'] + f'\n `[Source] <{prodinfo["link"]}>`_' )
            rst_output.append('parameters')
            params = parameters_to_dict(prodinfo['parameters'])
            for param, paraminfo in params.items():
                try:
                    rst_output.append(f'   {param.split()[0]} : *{param.split()[1]}*')
                except IndexError:
                    rst_output.append(f'   {param}')
                rst_output.append('       ' + paraminfo)
            rst_output.append('')
            # rst_output.append(f"\n\n{value2['parameters']}\n\n")
    return '\n'.join(rst_output)


def process_file(filepath):
    with open(filepath, 'r') as f:
        content = f.readlines()
    all_products = {}
    started = False
    product = None
    for i,line in enumerate(content):
        if 'def generate' in line:
            function = line.split('def')[1].split('(')[0].strip()
            all_products[function] = {}
            reading_title = False
            reading_params = False
        if 'Accepted product types:' in line:
            started = True
            reading_params = False
        if started:
            if 'All the products of the' in line:
                all_products_type = line.split('All the products of the ')[1].split()[0]
                all_products[function]['all_products'] = all_products_type

            match = re.findall("^'[A-Z0-9_]*'\\s*:\\d*", line.strip())
            if 'Parameters' in line: # End of block with product list
                reading_params = False
            if len(match):
                reading_params = False
                reading_title = True
                if product in all_products[function]:
                    all_products[function][product]['description'] = " ".join(
                        descr.replace('\n', '').split())

                product = match[0].replace("'", "").split(':')[0].strip()
                descr = line.split(':')[1].strip()
                all_products[function][product] = {}
                all_products[function][product]['parameters'] = ''
                continue
            if 'User defined parameters' in line:
                all_products[function][product]['description'] = " ".join(
                    descr.replace('\n', '').split())
                reading_params = True
                reading_title = False
                continue
            if reading_title:
                descr += line
            if reading_params and product:
                all_products[function][product]['parameters'] += " " + \
                    " ".join(line.replace('\n', ' ').split())
        if "prdcfg['type']" in line and '==' in line:
            for product in all_products[function].keys():
                if product in line:
                    all_products[function][product]['link'] = (funcpath_to_docpath(filepath) +
                        f'#L{i+1}')
    return all_products

products = {}
for root, _, files in os.walk(prodpath):
    for file in files:
        if file.endswith(".py") and file.startswith("process"):
            file_path = os.path.join(root, file)
            products.update(process_file(file_path))

# Match function names to product types
mapping_func_to_prodname = {}
file_path = Path(prodpath, 'product_aux.py')
with open(file_path, 'r') as f:
    content = f.readlines()
for line in content:
    match = re.findall("^'[A-Z0-9_]*'\\s*:\\d*", line.strip())
    if len(match):
        product = match[0].replace("'", "").split(':')[0]
        function = line.split(':')[1].strip()
        mapping_func_to_prodname[function] = product

for k in mapping_func_to_prodname:
    if k in products:
        products[mapping_func_to_prodname[k]] = products.pop(k)


# Convert dict data to reStructuredText format
rst_content = dict_to_restructured_text(products)
fname = Path(OUT_DIRECTORY, 'list_products.rst')
with open(fname, 'w') as f:
    f.write(rst_content)
