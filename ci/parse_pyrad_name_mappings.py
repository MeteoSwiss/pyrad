# %%
import pandas as pd
import inspect
from pathlib import Path

from pyrad.io import io_aux

FUNCTIONS_TO_PARSE = ['get_fieldname_pyart',
                      'get_datatype_odim',
                      'get_datatype_metranet',
                      'get_fieldname_icon']

mainpath = Path(__file__).resolve().parent.parent
OUT_DIRECTORY = str(Path(mainpath, 'doc', 'source', 'overview', 'mappings'))


for fun_name in FUNCTIONS_TO_PARSE:
    print('Parsing function {:s}'.format(fun_name))
    fun = getattr(io_aux, fun_name)
    mapping = fun_name.split('_')[-1]
    pyrad_dtypes = []
    srccode = inspect.getsource(fun)
    srccode = srccode.split('\n')
    for line in srccode:
        if ('datatype' in line or 'field_name' in line) and '==' in line:
            pyrad_dtypes.append(line.split('==')[1].split(':')[
                                0].strip().replace('"', '').replace("'",""))
        if 'return' in line:
            returnline = line.replace('return', '').strip()

    # Try to get output types from return statements of srccode
    all_values = []
    all_keys = []
    if '{' and '}' in returnline:
        keyname, valname = returnline.replace(
            '{', '').replace('}', '').split(':')
        keyname = keyname.strip()
        valname = valname.strip()
        if mapping not in keyname:
            keyname += '_' + mapping
    else:
        valname = returnline.strip()
    if mapping not in valname:
        valname += '_' + mapping

    for v in pyrad_dtypes:
        out = fun(v)
        if isinstance(out, dict):
            all_keys.append(list(out.keys())[0])
            all_values.append(list(out.values())[0])
        else:
            all_values.append(out)

    dic = {}
    dic['pyrad_name'] = pyrad_dtypes
    if len(all_keys):
        dic[keyname] = all_keys
    dic[valname] = all_values

    df = pd.DataFrame(dic)
    df.to_csv(str(Path(OUT_DIRECTORY, 'pyrad_to_{:s}.txt'.format(mapping))),
              index=False)
