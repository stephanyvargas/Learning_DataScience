from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import time

options = Options()
user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/33.0.1750.517 Safari/537.36'
options.add_argument('user-agent={0}'.format(user_agent))

driver = webdriver.Chrome(options=options, executable_path='/home/stephy/Selenium_driver_chrome/chromedriver_linux64/chromedriver')
wait = WebDriverWait(driver, 20)
action = ActionChains(driver)

driver.get("https://www.sigmaaldrich.com/JP/en/product/sigma/v5754")
html = driver.page_source
print(html)
#Login_Btn = wait.until(EC.element_to_be_clickable((By.XPATH, "//*[@class='pxc-fn-login']/a")))

#action.move_to_element(Login_Btn).click().perform()
