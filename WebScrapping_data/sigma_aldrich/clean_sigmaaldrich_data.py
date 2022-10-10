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


def get_unique_code(lst, how_many):
    code =[]
    for i in range(how_many):
        code.append(list(lst[i].keys())[0])
    return code


def get_available_products(lst, how_many, code):
    products = []
    for i in range(how_many):
        try:
            for product in lst[i][code[i]]['Available_products']:
                #j=i
                quantity = re.findall(r'\d+',product['ID - quatity'])[-1]
                if quantity and len(quantity) <= 5:
                    amount = quantity
                    unit = product['ID - quatity'].split(amount)[-1].strip()
                else:
                    amount, unit = None, None

                price = int(product['price'].replace(',',''))
                delivery_time = product['shipment']

                if 'Orders outside of US' or '米国および欧州外の注文の場合' in delivery_time:
                    stock_japan = False
                else:
                    stock_japan = True

                products.append({'code' : code[i],
                                 'amount' : amount,
                                 'unit' : unit,
                                 'price' : price,
                                 'stock_japan' : stock_japan,
                                 'aprox. delivery_time' : delivery_time})
        except ValueError:
            #print(f'WARNING ValueError: Check product {code[j]}, ', product)
            pass
        except TypeError:
            #print(f'WARNING TypeError: Check product {code[j]}, ', lst[j][code[j]])
            pass
        except IndexError:
            #print(f'WARNING IndexError: Check product {code[j]}, ', product)
            pass
        except KeyError:
            #print(f'WARNING KeyError: Check product {code[j]}, ', product)
            pass

    products_df = pd.DataFrame(products)
    products_df.index = products_df['code']
    products_df.drop(['code'], axis=1, inplace=True)
    return products_df


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

    num = len(entry)
    code = get_unique_code(entry, num)
    df_prod = get_available_products(entry, num, code)
    print(df_prod.sample(20))
    '''df = get_smiles(entry,num,code)
    print(df.info())
    print('original dataset : ', len(entry),
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
