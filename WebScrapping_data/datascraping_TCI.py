import requests
from bs4 import BeautifulSoup


response = requests.get(f"https://www.tcichemicals.com/JP/en/p/{A0001}", headers = {"User-Agent"
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


#Get the amount of product avaliable for purchase and price
size  = soup.find_all("td", attrs={"data-attr": "Size:"})
price = soup.find_all("td", attrs={"data-attr": "Unit Price"})

for a, p in zip(amount, price):
    a = re.sub('\\[a-z]','', a.string).strip()
    p = str(price[0]).split('¥')[-1].split('<')[0]
    info_dict[a] = p

