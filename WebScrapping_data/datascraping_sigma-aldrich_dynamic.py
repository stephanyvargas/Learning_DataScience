from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.proxy import *
from tqdm import tqdm
import datetime
import time
import json
import math
import yaml #print readeable dictionary in the terminal
from random import randint
import random


def colect_data(product_url, proxy_entry, product_data={}, inspect_page=False):
    #print(product_url,end=' ')

    #try:

    '''To not get detected -> and eventually blacklisted!
    # disabling enable-automation, or disabling automation controller disables
       webdriver.navigator which websites uses to detect automation scripts
    # Rotating the user-agent through execute_cdp_cmd()
    # Change the property value of the navigator for webdriver to undefined
    # Exclude the collection of enable-automation switches
    # Turn-off useAutomationExtension'''


    myProxy = 'http://{0}:{1}'.format(proxy_entry['IP Address'], proxy_entry['Port'])
    proxy = Proxy({
    'proxyType': ProxyType.MANUAL,
    'httpProxy': myProxy,
    'sslProxy': myProxy,
    'noProxy': ''})

    options = Options()
    #options.proxy = proxy
    options.add_argument('--headless')
    options.add_argument("start-maximized")
    options.add_argument('--no-sandbox')
    options.add_argument('--disable-setuid-sandbox')
    options.add_argument('--disable-dev-shm-usage')
    options.add_argument('--window-size=600,400')
    options.add_argument('--ignore-certificate-errors')
    options.add_argument('--disable-accelerated-2d-canvas')
    options.add_argument('--disable-gpu')
    options.add_argument('--proxy-server={}'.format(myProxy))
    options.add_experimental_option("excludeSwitches", ["enable-automation"])
    options.add_experimental_option('useAutomationExtension', False)
    driver = webdriver.Chrome(options=options,
                              executable_path='/home/stephy/Selenium_driver_chrome/chromedriver_linux64/chromedriver')
    driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")
    driver.execute_cdp_cmd('Network.setUserAgentOverride', {"userAgent": 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/83.0.4103.53 Safari/537.36'})


    ##executable_path is deprecated, need to change it eventually
    #Load the Chrome driver and open Chromium in the back.
    #A chrome tab will open and wait for the JavaScript to load.
    wait = WebDriverWait(driver, 40)
    action = ActionChains(driver)
    driver.get(product_url)


    #If the page source needs to be inspected
    if inspect_page:
        html = driver.page_source
        print(html)


    #####Load the embedded json in the webpage
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


    #Get the price, amount and shipping information.
    '''Use list_products instead of product_json["materialIds"] since there
    could be more products but are not listed as available'''
    try:
        availability = driver.find_element(By.CLASS_NAME, "MuiTableBody-root")
        list_products = availability.text.split('\n')
        product_availability = {}
        values = [val for val in list_products if val != '詳細...']
        num = len(values)/3
        var_a = 0
        CurrentDate = datetime.datetime.now()

        #Check that the shipping information table is properly built
        if math.ceil(num) != math.floor(num):
            #print('Number of entries is wrong!', ' Number:', num, math.ceil(num), math.floor(num))
            print('Number of entries is wrong! Check values for ', product_url)
            return False

        for product in range(int(num)):
            product_details = {}
            prod_var = values[var_a+0].split(' ')[0]
            product_details['ID - quatity'] = values[var_a+0]
            product_details['price'] = values[var_a+2].replace('￥', '')
            try:
                #If there is concrete information about the delivery date, get how many days it will take to deliver
                ExpectedDate = datetime.datetime.strptime(values[var_a+1].split(' ')[1], "%Y年%m月%d日")
                timedelta = ExpectedDate-CurrentDate
                product_details['shipment'] = str(abs(timedelta.days)) + ' day(s)'
            except:
                #If there is a message for the delivery instead of a date, save the message
                product_details['shipment'] = values[var_a+1]
            product_availability[prod_var] = product_details
            var_a += 3
        product_data['Available_products'] = product_availability
        product_availability = {}

    except:
        try:
            #In case there is no shipping information but there is a note
            elements = [element.get_attribute('innerHTML')\
                       for element in driver.find_element(By.ID, '__next').find_elements(By.TAG_NAME, "span")\
                       if 'Note' in element.get_attribute('innerHTML')]
            product_data['Available_products'] = {'Note' : elements[0].replace('<b>', '').replace('</b>', '').replace('Note: ', '') + 'Technical Service'}
        except:
            try:
                #In case the product is no longer available
                element = driver.find_element(By.TAG_NAME, 'strong').get_attribute('innerHTML')
                product_data['Available_products'] = {'WARNING' : 'Product {} might be discontinued.\
                                                       Contact Technical Service.'.format(element)}
            except:
                #When hope is lost!
                print('Something went wrong with innerHTML! Check {}'.format(product_url))
                #print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))
                return False

    #Once data has been fully loaded and scraped, close the chrome tab automatically.
    driver.close()

    #####To check the data that is being returned uncomment the next line
    #print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))
    return product_data



def get_proxy():
    #The list consists of ~100 IP addresses
    f = open('proxy_list.json')
    proxy_data = json.load(f)
    f.close()
    return proxy_data


dictionary = {}
proxy_list=get_proxy()

##If one particular element needs to be inspected:
#colect_data('https://www.sigmaaldrich.com/JP/en/product/sigald/179124',
#             proxy_entry=proxy_list[random.randint(0,len(proxy_list)-1)],
#             inspect_page=False)

with open('sigmaaldrich_products_urls.txt', 'r') as url_f, open("data_sigmaaldrich.json", 'a+') as data_f:
    urls_file = url_f.readlines()
    for url_ in tqdm(urls_file[:5]):
        url = url_.replace('\n', '')
        data = colect_data(url,
                           proxy_entry=proxy_list[random.randint(0,len(proxy_list)-1)])
        dictionary[url.split('/')[-2] + '_' + url.split('/')[-1]] = data
        if data:
            entry = json.load(data_f)
            entry.append(dictionary)
            file.seek(0)
            json.dump(entry, data_f, sort_keys=True, indent=4)
        #else:
        #    print("Could not save the data {}".format(url))
        dictionary = {}
        time.sleep(random.randint(5,10))
