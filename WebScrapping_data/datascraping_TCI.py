import requests
import re
import time
import pandas as pd
import json
from tqdm import tqdm
from bs4 import BeautifulSoup

code = 'A0001'

def colect_data(code, all_info={}):
    print(code,end=' ')

    if code in all_info.keys():
            pass

    else:
        #connect to the server and make it look like a normal connection through the browser
        response = requests.get("https://www.tcichemicals.com/JP/en/p/{}".format(code), headers = {"User-Agent":
                                "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"})
        soup = BeautifulSoup(response.content, "html.parser")


        #Dictionary where all information for this one chemical will be saved
        info_dict = {}

        title = soup.find("title").string.split('|')[0].strip()

        #Check if the code exists
        if "Page Not Found" in title:
            all_info[code] = {'error':'noInfo'}
            print(":NoInfo")
        else:
            print(':OK')


            try:
                #Get the name, CAS and code
                CAS = title.split(' ')[-1]
                name = title.strip(CAS).strip()
                info_dict['name'] = name
                info_dict['CAS'] = CAS
                info_dict['code'] = code
            except:
                pass


            try:
                #Get the amount of product avaliable for purchase, price and region
                amount  = soup.find_all("td", attrs={"data-attr": "Size:"})
                price   = soup.find_all("td", attrs={"data-attr": "Unit Price"})
                saitama = soup.find_all("td", attrs={"data-attr": "Saitama (Kawaguchi)"})
                hyogo   = soup.find_all("td", attrs={"data-attr": "Hyogo (Amagasaki)"})
                others  = soup.find_all("td", attrs={"data-attr": "Stock in other WH"})

                for i, a  in enumerate(amount):
                    #purchase information available per amount
                    info_per_amount={}

                    a = re.findall(r'\d+', str(a))[0]
                    p = str(price[i]).split('¥')[-1].split('<')[0]
                    info_per_amount['price'] = p

                    if 'Contact Us' in str(saitama[i]):
                        info_per_amount['Saitama (Kawaguchi)']= 'Contact Company'
                    else:
                        info_per_amount['Saitama (Kawaguchi)']= re.findall(r'\d+', str(saitama[i]))[0]

                    if 'Contact Us' in str(hyogo[i]):
                        info_per_amount['Hyogo (Amagasaki)']= 'Contact Company'
                    else:
                        info_per_amount['Hyogo (Amagasaki)']= re.findall(r'\d+', str(hyogo[i]))[0]

                    if 'Contact Us' in str(others[i]):
                        info_per_amount['Stock in other WH']= 'Contact Company'
                    else:
                        info_per_amount['Stock in other WH']= re.findall(r'\d+', str(others[i]))[0]

                    info_dict[a+'G'] = info_per_amount
            except:
                pass


            #Get the properties of the product
            table = soup.find_all("td")
            try:
                val = False
                for i, prop in enumerate(table):
                    if (prop.string is not None):
                        if prop.string == 'Product Number':
                            val = True
                        if val and i%2 == 0:
                            try:
                                key = re.sub('\\[a-z].','', prop.string).strip()
                                value = re.sub('\\[a-z].','', table[i+1].string).strip()
                                info_dict[key] = value.replace("\n", "")
                            except TypeError:
                                print('Error: ', str(prop))

            except:
                pass

            #Extract Molecular Formula  and the Molecular Weight.
            try:
                for i, prop in enumerate(table):
                    if (prop.string is not None) and ('Molecular Formula / Molecular Weight' in prop.string):
                        temp = table[i+1].text.split('\n')
                        info_dict['Molecular Formula'] = temp[0]
                        info_dict['Molecular Weight'] = temp[2].replace('=', '').strip()
            except:
                pass


            #Extract information from the Chemical Substances Control Law.
            try:
                for i, prop in enumerate(table):
                    if (prop.string is not None) and ('Chemical Substance Law' in prop.string):
                        temp = table[i+1].text.split('\n')
                        info_dict['Chemical Substance Law_No']= re.sub('\\[a-z].','', temp[-2]).strip()
            except:
                pass


            #Extract the compound category list
            try:
                category_list = soup.find_all("div", attrs={"class": "subCategory"})

                for cat in category_list:

                    rootCat = cat.find("h5").text
                    tree_list = cat.find_all("span",class_="startPoint")

                    catList = []
                    for tree in tree_list:
                        catList.append([c.text for c in tree.find_all("a")])
                    info_dict[rootCat] = catList
            except:
                pass

            #Get the shipment information
            try:
                ship = soup.find_all('div', class_='col-md-12 col-xs-12')
                info_dict['Shipment Information'] = ship[0].text.replace('\t', '').split('\n')[1]
            except:
                pass

            all_info[code] = info_dict

        return all_info



#Generate all the possible combinations of codes（A0000〜Z9999）
prefix_list = ['A']#[chr(i) for i in range(65,91)]
code_list = ['0001','0002','0003']#[str(s).zfill(4) for s in range(0,10000)]


#Since the data will be huge, save it as a JSON file
for prefix in prefix_list:

    all_info = {}

    for code_ in tqdm(code_list):
        code = prefix + code_

        try:
            all_info = colect_data(code, all_info)
        except:
            print(code + ' timeOut')
            all_info[code] = {'error':'timeout'}

        #Leave a 3 second interval to avoid overloading the server (avoid being blacklisted!)
        time.sleep(3)

df = pd.DataFrame(all_info.values(),index=all_info.keys())
df.to_json('data_TCI.json', orient = 'split', compression = 'infer', index = 'true')
