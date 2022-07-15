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
import yaml #print readeable dictionary in the terminal


#Create a user that would inspect the webpage (needed to avoid <Access Denied>)
options = Options()
user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/33.0.1750.517 Safari/537.36'
options.add_argument('user-agent={0}'.format(user_agent))


#Load the Chrome driver and open Chromium in the back.
##executable_path is deprecated, need to change it eventually
driver = webdriver.Chrome(options=options, executable_path='/home/stephy/Selenium_driver_chrome/chromedriver_linux64/chromedriver')
wait = WebDriverWait(driver, 10)
action = ActionChains(driver)


#A chrome tab will open and wait for the entire JavaScript to load.
driver.get("https://www.sigmaaldrich.com/JP/en/product/sigma/v5754")
#html = driver.page_source #Uncomment if the page source needs to be inspected


#Load the embedded json in the webpage
table = driver.find_element(By.ID, '__NEXT_DATA__')
convertedDict = json.loads(table.get_attribute('innerHTML'))

product_data = {}
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
availability = driver.find_element(By.CLASS_NAME, "MuiTableBody-root")
list_products = availability.text.split('\n')
if list_products:
    product_availability = {}
    num = len(list_products)/4
    var_a = 0
    CurrentDate = datetime. datetime. now()

    for product in range(int(num)):
        product_details = {}
        product_details['quatity'] = list_products[var_a+0].split(' ')[1] +\
                                     list_products[var_a+0].split(' ')[2]
        product_details['price'] = list_products[var_a+3]
        try:
            ExpectedDate = datetime.datetime.strptime(list_products[var_a+1].split(' ')[1], "%Y年%m月%d日")
            timedelta = ExpectedDate-CurrentDate
            product_details['shipment'] = str(abs(timedelta.days)) + ' days'
        except:
            product_details['shipment'] = list_products[var_a+1]
        product_availability[list_products[var_a+0].split(' ')[0]] = product_details
        var_a += 4

product_data['Available_products'] = product_availability

#Once data has been fully loaded and scraped, close the chrome tab automatically.
driver.close()

#check the file on the terminal
#print(yaml.dump(product_data, allow_unicode=True, default_flow_style=False))
return product_data
