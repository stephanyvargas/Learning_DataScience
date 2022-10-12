import pandas as pd
from rdkit import Chem
from utils import get_available_products, get_unique_code
import pprint
import json


def get_cansmi(smi):
    '''Input a Smiles String and tranform it to an RDKit'''
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
    except:
        return 'Error'


def get_attribute_list(lst,code):
    attr_list = []
    for i, compound in enumerate(lst):
        for attribute in lst[i][code[i]]['attributes']:
            if attribute['key'] not in attr_list:
                attr_list.append(attribute['key'])
            else: pass
    return attr_list


def get_smiles(lst,code):
    smiles = []
    for i, compound in enumerate(lst):
        for attribute in lst[i][code[i]]['attributes']:
            if attribute['key'] == 'smiles string':
                smile_string = attribute['values'][0]
                smiles.append({'code' : code[i],
                               'sigma_aldrich_smiles' : smile_string,
                               'rdkit_smiles' : None})

    smiles_df = pd.DataFrame(smiles)
    smiles_df['rdkit_smiles'] = smiles_df.sigma_aldrich_smiles.apply(get_cansmi)
    smiles_df.drop(['sigma_aldrich_smiles'], axis=1, inplace=True)
    error_df = smiles_df[smiles_df.rdkit_smiles == 'Error']
    return smiles_df


def print_table(df, sample=None, check_column=None):
    if sample:
        print(df.sample(20))
        print(df.info())
    elif check_column:
        print(df[check_column].unique())
    else:
        print(df)
        print(df.info())

def delete_unavailable_products(code, lst):
    for i, _ in enumerate(code):
        if i in lst:
            del code[i]
    return code


def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    products_list = entry#[:1000]


    # Extract the unique sigma aldrich identifier for each product
    code = get_unique_code(products_list)


    # Build availability table
    df_prod, df_naproduct = get_available_products(products_list, code)
    #print_table(df_prod, sample=20, check_column='special_remarks')
    #print(f'''scraped data: {len(code)}, number of unique codes: {len(set(code))}, database unique codes: {len(df_prod.code.unique())}''')


    ## Get rid of unavailable products, update the unique codes to only those that are available
    code = delete_unavailable_products(code, df_naproduct.idx)
    products_list = delete_unavailable_products(products_list, df_naproduct.idx)


    # Get the smiles for RDkit conversion
    print('smiles next')
    df_smiles = get_smiles(products_list,code)
    print_table(df_smiles)
    print('Products list: ', len(products_list))
    print('Smiles from sigma aldrich: ', len(df_smiles))
    print('Missing smiles representation: ', len(products_list) - len(df_smiles))


    # get available attributes
    #get_attribute_list(products_list,code)

if __name__ == "__main__":
    main()
