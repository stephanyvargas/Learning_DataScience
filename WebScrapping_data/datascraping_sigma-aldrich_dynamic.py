from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.proxy import *
from selenium.common.exceptions import WebDriverException
from tqdm import tqdm
import datetime
import time
import json
import math
import yaml #print readeable dictionary in the terminal
from random import randint
import random
import os


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
    options.proxy = proxy
    options.add_argument('--headless')
    options.add_argument("start-maximized")
    options.add_argument('--no-sandbox')
    options.add_argument('--disable-setuid-sandbox')
    options.add_argument('--disable-dev-shm-usage')
    options.add_argument('--window-size=600,400')
    options.add_argument('--ignore-certificate-errors')
    options.add_argument('--disable-accelerated-2d-canvas')
    options.add_argument('--disable-gpu')
    #options.add_argument('--proxy-server={}'.format(myProxy))
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

    try:
        driver.get(product_url)
    except WebDriverException:
        #skip if the internet connection is bad and is returning conectivity error
        print('Timeout')
        return False

    ####If the page source needs to be inspected
    if inspect_page:
        html = driver.page_source
        print(html)


    ####Load the embedded json in the webpage
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


    ####Get the price, amount and shipping information.
    '''Use list_products instead of product_json["materialIds"] since there
    could be more products but are not listed as available'''

    availability = False
    side_note = False
    no_longer_available = False
    not_available = False


    ####Load the shipping information table (if it exists)
    try:
        availability = driver.find_element(By.CLASS_NAME, "MuiTableBody-root")
    except:
        pass

    ####In case there is/isn't shipping information and/but there is a note
    try:
        driver.execute_script("window.scrollTo(0, 100)")
        side_note = [element.get_attribute('innerHTML')\
                    for element in driver.find_element(By.ID, '__next').find_elements(By.TAG_NAME, "span")\
                    if 'Note' in element.get_attribute('innerHTML')]
    except:
        pass

    ####In case the product is no longer available
    try:
        no_longer_available = driver.find_element(By.TAG_NAME, 'strong').get_attribute('innerHTML')
    except:
        pass

    ####In case the product is currently unavailable or is not for sale in Japan
    try:
        not_available = driver.find_element(By.XPATH, '//*[@id="prodductDetailGrid"]/div[1]/div/div/div[2]/div/div[3]/div/div/div/div/div/span[1]').get_attribute('innerHTML')
    except:
        pass


    if availability:
        list_products = availability.text.split('\n')
        product_availability = []
        values = [val for val in list_products if val != '詳細...']
        idx = [i for i, value in enumerate(values) if values[0][:4] in value]
        idx_values = [idx[i] - idx[i-1] for i in range(1, len(idx))]
        idx_values.append(len(values)-idx[-1])
        var_a = 0
        for n, product in enumerate(idx):
            product_details = {}
            prod_var = values[var_a+0].split(' ')[0]
            product_details['ID - quatity'] = values[var_a+0]
            product_details['shipment'] = get_delivery_time(values, var_a)
            if idx_values[n] == 3:
                product_details['price'] = values[var_a+2].replace('￥', '')
            if idx_values[n] == 4:
                product_details['price'] = values[var_a+3].replace('￥', '')
                product_details['note'] = values[var_a+2].replace('￥', '')
            product_availability.append(product_details)
            var_a += idx_values[n]
        product_data['Available_products'] = product_availability

    elif side_note:
        product_data['Side Note'] = {'Note':side_note[0].replace('<b>', '').replace('</b>', '').replace('Note: ', '') + 'Technical Service'}
        if not availability:
            product_data['Available_products'] = None

    elif no_longer_available:
        product_data['Side Note'] = {'WARNING' : 'Product {} might be discontinued. Contact Technical Service.'.format(no_longer_available)}
        if not availability:
            product_data['Available_products'] = None

    elif not_available:
        product_data['Side Note'] = {'WARNING' : not_available}
        if not availability:
            product_data['Available_products'] = None

    else:
        print('Product info is not yet supported by this code', product_url)
        return False

    #####To check the data that is being returned uncomment the next line
    #print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))
    driver.close()
    return product_data


def get_delivery_time(values, var_a):
    try:
        #If there is concrete information about the delivery date, get how many days it will take to deliver
        CurrentDate = datetime.datetime.now()
        ExpectedDate = datetime.datetime.strptime(values[var_a+1].split(' ')[1], "%Y年%m月%d日")
        timedelta = ExpectedDate-CurrentDate
        return str(abs(timedelta.days)) + ' day(s)'
    except:
        #If there is a message for the delivery instead of a date, save the message
        return values[var_a+1]


def get_proxy():
    #The list consists of ~100 IP addresses, can be updated running the get_proxy.py script
    f = open('proxy_list.json')
    proxy_data = json.load(f)
    f.close()
    return proxy_data


dictionary = {}
proxy_list=get_proxy()

##If one particular element needs to be inspected:
#colect_data('https://www.sigmaaldrich.com/JP/en/product/aldrich/eme00109',
#             proxy_entry=proxy_list[random.randint(0,len(proxy_list)-1)],
#             inspect_page=False)

with open('sigmaaldrich_products_urls.txt', 'r') as url_f:
    urls_file = url_f.readlines()
    for url_ in tqdm(urls_file[:2000]):
        url = url_.replace('\n', '')
        data = colect_data(url,
                           proxy_entry=proxy_list[random.randint(0,len(proxy_list)-1)])
        dictionary[url.split('/')[-2] + '_' + url.split('/')[-1]] = data
        if data:
            with open("data_sigmaaldrich.json", 'r+', encoding='utf8') as data_f:
                entry = json.load(data_f)
                entry.append(dictionary)
                data_f.seek(0)
                json.dump(entry, data_f, indent=4)
            entry.clear()
        dictionary = {}
        time.sleep(random.randint(5,10))
