import pandas as pd
from rdkit import Chem
import pprint
import json
import re

def get_cansmi(smi):
        try:
            return Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
        except:
            return None

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


def get_available_products(lst, code):
    products = []
    na_products = []
    for i, compound in enumerate(lst):
        try:
            # Get warnings from the company regarding products
            if hasattr(compound[code[i]]['Available_products'], 'get'):
                warning = compound[code[i]]['Available_products'].get('WARNING')
                title = compound[code[i]]['Metadata'].get('title')
                url = compound[code[i]]['url']
                na_products.append({'code' : code[i], 'title' : title, 'url' : url, 'warning' : warning})
                pass
            elif compound[code[i]]['Available_products'] == None:
                warning = compound[code[i]].get('Side Note')
                na_products.append({'code' : code[i], 'title' : None, 'url' : None, 'warning' : warning})
                pass
            else:
                #print(compound[code[i]])
                for product in compound[code[i]]['Available_products']:

                    # Check if there are any special notes on the product
                    note = product.get('note')

                    # Get the quatity and the units of the product.
                    quantity = re.findall(r'\d+',product['ID - quatity'])[-1]
                    if quantity and len(quantity) <= 5:
                        amount = quantity
                        unit = product['ID - quatity'].split(amount)[-1].strip()
                    else:
                        amount, unit = None, None

                    # Check if there are products that need to be shipped from USA
                    delivery_time = product['shipment']
                    if '申し訳ございませんが、この製品のフルフィルメントと配送が遅延しています。' in delivery_time:
                        delivery_time = 'Delayed'
                    if ('Orders outside of US' in delivery_time) or ('米国および欧州外の注文の場合' in delivery_time):
                        stock_japan = False
                    else:
                        stock_japan = True

                    # Price might not be available and the company might need to be contacted.
                    if product['price'] == 'お問い合わせ':
                        price = None
                    else:
                        price = int(product['price'].replace(',',''))

                    products.append({'code' : code[i],
                                     'amount' : amount,
                                     'unit' : unit,
                                     'price' : price,
                                     'stock_japan' : stock_japan,
                                     'aprox_delivery_time' : delivery_time,
                                     'special_remarks' : note})

        except ValueError:
            print(f'WARNING ValueError: Check compound {code[i]}, ', compound)
            pass
        except TypeError:
            print(f'WARNING TypeError: Check compound {code[i]}, ', compound)
            pass
        except IndexError:
            print(f'WARNING IndexError: Check compound {code[i]}, ', compound)
            pass
        except KeyError:
            print(f'WARNING KeyError: Check compound {code[i]}, ', compound)
            pass

    products_df = pd.DataFrame(products)
    return products_df, na_products


def get_smiles(lst,how_many,code):
    smiles = []
    for i in range(how_many):
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
    return smiles_df

def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    products_list = entry[1850:1880]
    code = get_unique_code(products_list)
    df_prod, discont_lst = get_available_products(products_list, code)
    print(df_prod)
    print(discont_lst)
    '''df = get_smiles(products_list,num,code)
    print(df.info())
    print('original dataset : ', len(products_list),
          'available products : ', len(df_prod),
          'smiles dataset : ', len(df))
    '''
    #print(len(entry))

if __name__ == "__main__":
    main()






#get_duplicates(lst)
##if there are duplicates, eliminate through:
#with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
#    entry = json.load(data_f)
#    del entry[3] ##change IT to n!!!!!!
#    data_f.seek(0)
#    json.dump(entry, data_f, sort_keys=True, indent=4)
