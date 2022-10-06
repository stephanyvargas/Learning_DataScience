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


def main():
    with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
        entry = json.load(data_f)

    print(get_unique_code(entry, 20))

    '''    for product in entry[i][code[i]]['Available_products']:
            quantity = re.findall(r'\d+',product['ID - quatity'])[-1]
            units = product['ID - quatity'].split(quantity)[-1].strip()
            print(i, ': ', quantity, units)
    '''
    #print(code)

if __name__ == "__main__":
    main()






#get_duplicates(lst)
##if there are duplicates, eliminate through:
#with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
#    entry = json.load(data_f)
#    del entry[3] ##change IT to n!!!!!!
#    data_f.seek(0)
#    json.dump(entry, data_f, sort_keys=True, indent=4)

'''
with open('sigmaaldrich_products_urls.txt', 'r') as url_f:
    urls_file = url_f.readlines()

with open('sigmaaldrich_products_ALLurls.txt', 'r') as ALLurl_f:
    ALLurls_file = ALLurl_f.readlines()

indermediate_urls = [u.replace('\n', '') for u in urls_file]
not_done = [prod for prod in indermediate_urls if prod not in lst]

print('not done:', len(not_done),'done:', len(lst), 'indermediate_urls:', len(indermediate_urls))
print('All products = Not Done + Done:', len(ALLurls_file), '=', len(not_done) + len(lst), '-> these numbers should be the same!')

with open('sigmaaldrich_products_urls_intermediate.txt', 'w') as fp:
    for item in not_done:
        # write each item on a new line
        fp.write("%s\n" % item)
    print('Done')
'''
