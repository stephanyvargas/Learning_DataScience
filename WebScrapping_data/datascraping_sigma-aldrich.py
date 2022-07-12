import requests
import re
import time
import pandas as pd
import json
from tqdm import tqdm
from bs4 import BeautifulSoup
import os

urls_file = open("sigmaaldrich_products_urls.txt", "r")
urls_list = urls_file.read().split('\n')

one_url=urls_list[0]
print(one_url)

response = requests.get(one_url, headers = {"User-Agent":
                        "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:47.0) Gecko/20100101 Firefox/47.0"})
soup = BeautifulSoup(response.content, "html.parser")
print(soup)
