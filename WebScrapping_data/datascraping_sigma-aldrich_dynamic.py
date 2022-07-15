from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import time


#Create a user that would inspect the webpage (needed to avoid <Access Denied>)
options = Options()
user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/33.0.1750.517 Safari/537.36'
options.add_argument('user-agent={0}'.format(user_agent))

#Load the Chrome driver and open Chromium in the back.
driver = webdriver.Chrome(options=options, executable_path='/home/stephy/Selenium_driver_chrome/chromedriver_linux64/chromedriver')
wait = WebDriverWait(driver, 10)
action = ActionChains(driver)

#A chrome tab will open and wait for the entire JavaScript to load.
driver.get("https://www.sigmaaldrich.com/JP/en/product/sigma/v5754")
html = driver.page_source
heading1 = driver.find_element(By.CLASS_NAME, "MuiTableBody-root")
print(heading1.text)
#print(html)

#Once data has been fully loaded and scraped, close the chrome tab automatically.
driver.close()
