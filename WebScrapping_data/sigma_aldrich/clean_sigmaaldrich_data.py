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
                url = compound[code[i]].get('url')
                na_products.append({'code' : code[i], 'title' : title, 'url' : url, 'warning' : warning})

            elif compound[code[i]]['Available_products'] == None:
                side_note = compound[code[i]].get('Side Note')
                note = side_note.get('WARNING')
                maybe_available = side_note.get('Note')

                if maybe_available:
                    # This product can be packaged on demand. Need to contact company directly.
                    products.append({'index' : i,
                                     'code' : code[i],
                                     'amount' : None,
                                     'unit' : None,
                                     'price' : None,
                                     'stock_japan' : None,
                                     'aprox_delivery_time' : None,
                                     'special_remarks' : maybe_available})


                elif ('現在、価格および在庫状況を閲覧できません' in note) or ('Pricing and availability is not currently available' in note):
                    # Need to contact company directly for information.
                    products.append({'index' : i,
                                     'code' : code[i],
                                     'amount' : None,
                                     'unit' : None,
                                     'price' : None,
                                     'stock_japan' : None,
                                     'aprox_delivery_time' : None,
                                     'special_remarks' : note})
                    

                title = compound[code[i]]['Metadata'].get('title')
                url = compound[code[i]].get('url')
                na_products.append({'code' : code[i], 'title' : title, 'url' : url, 'note' : note})

            else:
                for j, product in enumerate(compound[code[i]]['Available_products']):

                    # Check if there are any special notes on the product
                    note = product.get('note')
                    if note and (note[-1] == '0'):
                        note == None

                    # Get the quatity and the units of the product.
                    quantity = re.findall(r'\d+',product['ID - quatity'])
                    if quantity and len(quantity) <= 5:
                        amount = quantity[-1]
                        unit = product['ID - quatity'].split(amount)[-1].strip()
                    else:
                        amount, unit = None, None

                    # Check if there are products that need to be shipped from USA
                    delivery_time = product.get('shipment')
                    if '申し訳ございませんが、この製品のフルフィルメントと配送が遅延しています。' in delivery_time:
                        delivery_time = 'Delayed'
                    elif ('Orders outside of US' in delivery_time) or ('米国および欧州外の注文の場合' in delivery_time):
                        stock_japan = False
                        delivery_time = 'Delivered from the US, may take several weeks.'
                    elif delivery_time:
                        stock_japan = True
                    else:
                        stock_japan = False

                    # Price might not be available and the company might need to be contacted.
                    if (product.get('price') == 'お問い合わせ') or not product.get('price'):
                        price = None
                    else:
                        price = int(product['price'].replace(',',''))

                    # Error from the scraping algorithm, same product, same units and amount mixed.
                    if (delivery_time[-1] == '0') or (',' in delivery_time) and not products[-1].get('price'):
                        print('last product: ', products[-1])
                        print('current product: ', product)
                        products[-1]['price'] = int(delivery_time.replace(',','')[1:])
                        print('last -1 price: ', products[-1]['price'])# = int(delivery_time.replace(',',''))

                    else:
                        products.append({'index' : i,
                                         'code' : code[i],
                                         'amount' : amount,
                                         'unit' : unit,
                                         'price' : price,
                                         'stock_japan' : stock_japan,
                                         'aprox_delivery_time' : delivery_time,
                                         'special_remarks' : note})

        except ValueError:
            print(i, code[i])
            print(f'WARNING ValueError: Check compound {code[i]}, ', compound)
            pass
        except TypeError:
            print(i, code[i])
            print(f'WARNING TypeError: Check compound {code[i]}, ', compound)
            pass
        except IndexError:
            print(i, code[i])
            print(f'WARNING IndexError: Check compound {code[i]}, ', compound)
            pass
        except KeyError:
            print(i, code[i])
            print(f'WARNING KeyError: Check compound {code[i]}, ', compound)
            pass

    products_df = pd.DataFrame(products)
    na_products_df = pd.DataFrame(na_products)
    return products_df, na_products_df


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

    products_list = entry#[12070:12080]
    code = get_unique_code(products_list)
    #get_available_products(products_list, code)
    df_prod, df_naproduct = get_available_products(products_list, code)
    print(df_prod)
    #print(df_naproduct.note.unique())
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
