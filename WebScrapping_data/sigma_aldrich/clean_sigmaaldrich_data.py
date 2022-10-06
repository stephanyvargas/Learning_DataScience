import pandas as pd
import json
import re

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
        for product in lst[i][code[i]]['Available_products']:
            try:
                quantity = re.findall(r'\d+',product['ID - quatity'])[-1]
                if quantity and len(quantity) <= 5:
                    amount = quantity
                    unit = product['ID - quatity'].split(amount)[-1].strip()
                else:
                    amount, unit = None, None

                price = int(product['price'].replace(',',''))
                delivery_time = product['shipment']

                if 'Orders outside of US' in delivery_time:
                    stock_japan = False
                else:
                    stock_japan = True

                products.append({'code' : code[i],
                                 'amount' : amount,
                                 'unit' : unit,
                                 'price' : price,
                                 'stock_japan' : stock_japan,
                                 'aprox. delivery_time' : delivery_time})
            except TypeError:
                print(f'WARNING: Check product {code[i]}')
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
                smiles.append({'code' : code[i],
                               'sigma_aldrich_smiles' : smile_string})
    smiles_df = pd.DataFrame(smiles)
    smiles_df.index = smiles_df['code']
    smiles_df.drop(['code'], axis=1, inplace=True)
    return smiles_df


def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    num = 10
    code = get_unique_code(entry, num)
    #print(get_available_products(entry, num, code))
    print(get_smiles(entry,num,code))

if __name__ == "__main__":
    main()






#get_duplicates(lst)
##if there are duplicates, eliminate through:
#with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
#    entry = json.load(data_f)
#    del entry[3] ##change IT to n!!!!!!
#    data_f.seek(0)
#    json.dump(entry, data_f, sort_keys=True, indent=4)
