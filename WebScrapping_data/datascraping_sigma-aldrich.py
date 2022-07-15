import json
import requests
from bs4 import BeautifulSoup
#import yaml #print pretty dictionary in the terminal
#print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))


'''Warning! This version does not include price or date availability!!'''


def colect_data(url_prod):

    #Connection to the server
    webpage = requests.get(url_prod, headers = {"User-Agent":
                            "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"})
    soup = BeautifulSoup(webpage.content, "html.parser")

    #Initialize an empty dictionary to save the data
    product_data = {}

    try:
        #Get the available Json information and use as a dictionary
        table = soup.find(id = '__NEXT_DATA__')
        convertedDict = json.loads(table.text)
        #print(yaml.dump(convertedDict, allow_unicode=True, default_flow_style=False))

        product_json = convertedDict["props"]["pageProps"]["data"]["getProductDetail"]
        product_data['id'] = product_json['id']
        product_data['brand'] = product_json["brand"]["name"]
        product_data['aliases'] = product_json["aliases"]
        product_data['name'] = product_json["name"]
        product_data["molecularWeight"] = product_json["molecularWeight"]
        product_data["empiricalFormula"] = product_json["empiricalFormula"]
        product_data["linearFormula"] = product_json["linearFormula"]
        product_data["synonyms"] = product_json["synonyms"]
        product_data["attributes"] = product_json["attributes"] #Here is the SMILES string!
        product_data["materialIds"] = product_json["materialIds"] #Available products, no amount, price or shipping info available...
        product_data["compliance"] = product_json["compliance"]
        product_data["complianceReach"] = product_json["complianceReach"]
        product_data["Metadata"] = product_json["browserMetadata"] #Basic information displayed in the website
        product_data["features"] = product_json["features"]
        product_data["components"] = product_json["components"]
        product_data["substanceCount"] = product_json["substanceCount"]
        product_data["productCategories"] = product_json["productCategories"]
        product_data["catalogId"] = product_json["catalogId"]
        product_data["paMessage"] = product_json["paMessage"]
        product_data["relatedProducts"] = product_json["relatedProducts"]
        product_data["type"] = product_json["type"]
        return product_data['id'], product_data

    except:
        print('Something went wrong here! Check {}'.format(url_prod))
        return False


def available_products_url():
    #The current list consists of ~100_000 products as of July 2022
    ##Open the list of available products urls scraped using the script "urlscraping_sigma-aldrich.py"
    urls_file = open("sigmaaldrich_products_urls.txt", "r")
    urls_list = urls_file.read().split('\n')
    return urls_list[:2] #while testing return few products


#Get all the possible products from the system.
list_products = available_products_url()


#Get the individual information for each product
for product in list_products:

    all_info = {}
    try:
        product_id, product_information = colect_data(product)
        all_info[product_id] = product_information
    except:
        print(product + ' timeOut')

print(all_info)

##Save the data
#with open('data.json', 'w') as fp:
#    json.dump(all_info, fp,  indent=4)

    #Leave a 3 second interval to avoid overloading the server (avoid being blacklisted!)
time.sleep(3)
#
#    if all_info:
#        path = 'data_sigma_aldrich.json'
#        try:
#            with open('data.json', 'w') as fp:
#                json.dump(all_info, fp,  indent=4)
#
#        except:
#            print('Data for codes list {} could not be saved'.format(prefix))
#
