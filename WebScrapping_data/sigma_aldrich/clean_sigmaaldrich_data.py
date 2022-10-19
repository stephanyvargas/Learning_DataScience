import pandas as pd
from utils import get_available_products, get_unique_code
from utils import get_smiles, save_smiles_txt
from utils import delete_unavailable_products
import json
import re


def get_aliases_list(lst,code,save_output=False,name='name',directory=None):
    alias_df = pd.DataFrame()
    for i, compound in enumerate(lst):
        for attribute in lst[i][code[i]]['aliases']:
            if attribute['key'] == 'pubchem substance id':
                rep = r'>(\d+)<'
                value = re.findall(rep, attribute['value'])[0]
                alias_df.loc[code[i], attribute['key']]=value
            else:
                alias_df.loc[code[i], attribute['key']]=attribute['value']
    
    if save_output:
        save_df_json(alias_df, name)
    return alias_df


def save_df_json(df, name, directory=None):
    if directory == None:
        directory = ''
    if directory and directory[-1] != '/':
        directory += '/'
    else: pass
    path = f'{directory}{name}.json'
    df.to_json(path, orient = 'split', compression = 'infer', index = 'true')


def print_table(df, sample=None, check_column=None):
    if sample:
        print(df.sample(20))
        print(df.info())
    elif check_column:
        print(df[check_column].unique())
    else:
        print(df)
        print(df.info())


def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    products_list = entry#[:50]
    scraped_products = len(products_list)

    # Extract the unique sigma aldrich identifier for each product
    code = get_unique_code(products_list)


    '''Build availability table'''
    df_prod, df_naproduct = get_available_products(products_list, code)
    #save_df_json(df_prod, 'sigma_aldrich_package')
    #print_table(df_prod, sample=20, check_column='special_remarks')


    '''Get rid of unavailable products, update the unique codes to only those that are available'''
    code = delete_unavailable_products(code, df_naproduct.idx)
    products_list = delete_unavailable_products(products_list, df_naproduct.idx)


    '''Get the smiles for RDkit conversion'''
    #df_smiles, df_smiles_error = get_smiles(products_list,code,save_output=False, name='sigma_aldrich')
    #print_table(df_smiles, sample=20)


    '''Get the compound identifications'''
    print(get_aliases_list(products_list,code))


    '''Get the listed attributes'''
    #print(get_attribute_list(products_list,code))



    '''
    print('Scraped compounds list: ', scraped_products)
    print('Purchasable compounds list: ', len(products_list))
    print('Available products to purchase: ', len(df_prod))
    print('Smiles from sigma aldrich: ', len(df_smiles))
    print('Missing smiles representation: ', len(products_list) - len(df_smiles))
    '''


if __name__ == "__main__":
    main()
