import json

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


with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
    entry = json.load(data_f)

lst = []
for product in entry:
    prod = list(product.keys())[0].replace('_', '/')
    lst.append('https://www.sigmaaldrich.com/product/'+prod)

get_duplicates(lst)
