import pandas as pd
from rdkit import Chem
from utils import get_available_products
import pprint
import json
import re

def get_cansmi(smi):
    '''Input a Smiles String and tranform it to an RDKit'''
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
    except:
        print(smi)
        return None#'Error'

def get_duplicates(lst):
    #check for repeated elements
    newlist = []
    duplist = []
    for n, i in enumerate(lst):
        if i not in newlist:
            newlist.append(i)
        else:
            duplist.append(i)
            print('repeated element index:', n, 'url:', i)
            return duplist


def get_unique_code(lst):
    code =[]
    for value in lst:
        code.append(list(value.keys())[0])
    return code


def get_smiles(lst,code):
    smiles = []
    for i, compound in enumerate(lst):
        for attribute in lst[i][code[i]]['attributes']:
            if attribute['key'] == 'smiles string':
                smile_string = attribute['values'][0]
                smile_rdkit = None
                smiles.append({'code' : code[i],
                               'sigma_aldrich_smiles' : smile_string,
                               'rdkit_smiles' : smile_rdkit})
    smiles_df = pd.DataFrame(smiles)
    smiles_df.index = smiles_df['code']
    smiles_df['rdkit_smiles'] = smiles_df.sigma_aldrich_smiles.apply(get_cansmi)
    smiles_df.drop(['code'], axis=1, inplace=True)
    #print(smiles_df[smiles_df.rdkit_smiles == 'Error'])
    return smiles_df

def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    products_list = entry#[:1000]
    code = get_unique_code(products_list)

    # Build availability table
    #df_prod, df_naproduct = get_available_products(products_list, code)
    print_available_table = False
    if print_available_table:
        print(df_prod.sample(20))
        print(df_prod.info())
        print(df_naproduct)
        print(df_prod.special_remarks.unique())
        print(f'scraped data: {len(code)}, number of unique codes: {len(set(code))}, database unique codes: {len(df_prod.code.unique())}')

    # Get the smiles for RDkit conversion
    # --------> filter out the df_naproducts
    print('smiles next')
    df = get_smiles(products_list,code)
    print(df.info())
    '''print(df.info())
    print('original dataset : ', len(products_list),
          'available products : ', len(df_prod),
          'smiles dataset : ', len(df))
    '''
    #print(len(entry))

if __name__ == "__main__":
    main()
