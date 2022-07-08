import requests
import re
import time
from bs4 import BeautifulSoup

def collect_urls(page_num=1, page = True):
    products_list = []

    while page:
        search_url = 'https://www.sigmaaldrich-jp.com/structure/?page={}'.format(page_num)
        response = requests.get(search_url, headers = {"User-Agent":
                                    "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"})
        soup = BeautifulSoup(response.content, "html.parser")

        title = soup.find("title").string.split('|')[0].strip()

        if "Error" in title:
            #print(":NoInfo")
            page=False
            break
        else:
            print(':OK')

        products_url = soup.find_all('div', class_='products-name')

        for i, prod in enumerate(products_url):
            product = str(prod).split('"')[5]
            products_list.append(product)
            print(product)

        print(page_num)
        page_num += 1

        #Leave a 3 second interval to avoid overloading the server (avoid being blacklisted!)
        time.sleep(3)

    return products_list

products = collect_urls()#page_num=5108)

with open('sigmaaldrich_products_urls.txt', 'w') as fp:
    for item in products:
        # write each item on a new line
        fp.write("%s\n" % item)
    print('Done')
