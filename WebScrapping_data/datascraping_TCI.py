import requests
from bs4 import BeautifulSoup


response = requests.get("https://www.tcichemicals.com/JP/en/p/A0001", headers = {"User-Agent"
                        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"})

soup = BeautifulSoup(response.content, "html.parser")


title = soup.find("title").string.split('|')[0].strip()
CAS = title.split(' ')[-1]
name = title.strip(CAS).strip()


size  = soup.find_all("td", attrs={"data-attr": "Size:"})
price = soup.find_all("td", attrs={"data-attr": "Unit Price:"})
