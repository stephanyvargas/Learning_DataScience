import requests
import re
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


        #Get the name, CAS and code
        title = soup.find("title").string.split('|')[0].strip()
        CAS = title.split(' ')[-1]
        name = title.strip(CAS).strip()

        info_dict['name'] = name
        info_dict['CAS'] = CAS
        info_dict['code'] = code


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


        #Get the properties of the product
        table = table = soup.find_all("td")

        for i, prop in enumerate(table):
            if (prop.string is not None) and (('Melting point' in prop.string) or \
            ('Boiling point' in prop.string) or ('Purity / Analysis Method' in prop.string) or \
            ('PubChem Substance ID' in prop.string) or ('Merck Index (14)' in prop.string) or \
            ('Refractive index' in prop.string) or ('Density' in prop.string) or \
            ('Dissolution' in prop.string) or ('Flash point' in prop.string) or \
            ('Appearance' in prop.string) or ('Dissolution' in prop.string) or \
            ('Poisonous and Deleterious Substances' in prop.string) or ('Maximum Absorption' in prop.string) ):

                key = re.sub('\\[a-z].','', prop.string).strip()
                value = re.sub('\\[a-z].','', table[i+1].string).strip()
                info_dict[key] = value.replace("\n", "")
                
        return info_dict
