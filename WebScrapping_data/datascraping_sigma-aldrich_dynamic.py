from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import datetime
import time
import json
import math
import yaml #print readeable dictionary in the terminal


def colect_data(product_url, product_data={}):
    print(product_url,end=' ')

    #try:
    #Create a user that would inspect the webpage (needed to avoid <Access Denied>)
    options = Options()
    user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/33.0.1750.517 Safari/537.36'
    options.add_argument('user-agent={0}'.format(user_agent))


    #Load the Chrome driver and open Chromium in the back.
    ##executable_path is deprecated, need to change it eventually
    driver = webdriver.Chrome(options=options, executable_path='/home/stephy/Selenium_driver_chrome/chromedriver_linux64/chromedriver')
    wait = WebDriverWait(driver, 25)
    action = ActionChains(driver)


    #A chrome tab will open and wait for the entire JavaScript to load.
    driver.get(product_url)
    #html = driver.page_source #Uncomment if the page source needs to be inspected


    #Load the embedded json in the webpage
    table = driver.find_element(By.ID, '__NEXT_DATA__')
    convertedDict = json.loads(table.get_attribute('innerHTML'))

    product_json = convertedDict["props"]["pageProps"]["data"]["getProductDetail"]
    product_data["url"] = product_url
    product_data['id'] = product_json['id']
    product_data['brand'] = product_json["brand"]["name"]
    product_data['aliases'] = product_json["aliases"]
    product_data['name'] = product_json["name"]
    product_data["molecularWeight"] = product_json["molecularWeight"]
    product_data["empiricalFormula"] = product_json["empiricalFormula"]
    product_data["linearFormula"] = product_json["linearFormula"]
    product_data["synonyms"] = product_json["synonyms"]
    product_data["attributes"] = product_json["attributes"] #Here is the SMILES string!
    product_data["materialIds"] = product_json["materialIds"] #Available products, no amount, price or shipping info available.
    product_data["compliance"] = product_json["compliance"]
    product_data["complianceReach"] = product_json["complianceReach"]
    product_data["Metadata"] = product_json["browserMetadata"] #Basic information displayed in the website, Includes CAS if available
    product_data["features"] = product_json["features"]
    product_data["components"] = product_json["components"]
    product_data["substanceCount"] = product_json["substanceCount"]
    product_data["productCategories"] = product_json["productCategories"]
    product_data["catalogId"] = product_json["catalogId"]
    product_data["paMessage"] = product_json["paMessage"]
    product_data["relatedProducts"] = product_json["relatedProducts"]
    product_data["type"] = product_json["type"]


    '''Use list_products instead of product_json["materialIds"] since there
    could be more products but are not listed as available'''
    #Get the price, amount and shipping information.
    try:
        availability = driver.find_element(By.CLASS_NAME, "MuiTableBody-root")
        print('it is here!')
        list_products = availability.text.split('\n')
        product_availability = {}
        values = [val for val in list_products if val != '詳細...']
        num = len(values)/3
        var_a = 0
        CurrentDate = datetime.datetime.now()

        if math.ceil(num) != math.floor(num):
            print('Number of entries is wrong!', ' Number:', num, math.ceil(num), math.floor(num))

        for product in range(int(num)):
            product_details = {}
            prod_var = values[var_a+0].split(' ')[0]
            product_details['quatity'] = prod_var.rsplit('-', 1)[1]
            product_details['price'] = values[var_a+2]
            try:
                ExpectedDate = datetime.datetime.strptime(values[var_a+1].split(' ')[1], "%Y年%m月%d日")
                timedelta = ExpectedDate-CurrentDate
                product_details['shipment'] = str(abs(timedelta.days)) + ' day(s)'
            except:
                product_details['shipment'] = values[var_a+1]
            product_availability[prod_var] = product_details
            var_a += 3
        product_data['Available_products'] = product_availability

    except:
        try:
            elements = [element.get_attribute('innerHTML')\
                       for element in driver.find_element(By.ID, '__next').find_elements(By.TAG_NAME, "span")\
                       if 'Note' in element.get_attribute('innerHTML')]
            product_data['Available_products'] = elements[0].replace('<b>', '').replace('</b>', '')
        except:
            pass

    #Once data has been fully loaded and scraped, close the chrome tab automatically.
    driver.close()
    #print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))

    #except:
    #    print('Something went wrong here! Check {}'.format(product_url))
    #    print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))
    #    return False
    return product_data


dictionary = {}

colect_data('https://www.sigmaaldrich.com/JP/en/product/aldrich/798096')

#with open('sigmaaldrich_products_urls.txt', 'r') as url_f, open("data_sigmaaldrich.json", 'a+') as data_f:
#    urls_file = url_f.readlines()
#    for i, url in enumerate(urls_file[:100]):
#        data = colect_data(url)
#        dictionary[url.split('/')[-2] + '_' + url.split('/')[-1]] = data
#        print(i+1)
#        if data:
#            json.dump(dictionary, data_f, sort_keys=True, indent=4)
#
#        dictionary = {}
#
