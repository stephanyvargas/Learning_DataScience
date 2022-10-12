import pandas as pd
import re

def get_unique_code(lst):
    code =[]
    for value in lst:
        code.append(list(value.keys())[0])
    return code

def get_available_products(lst, code):
    products = []
    na_products = []
    for i, compound in enumerate(lst):
        try:
            # Get warnings from the company regarding products
            if hasattr(compound[code[i]]['Available_products'], 'get'):
                warning = compound[code[i]]['Available_products'].get('WARNING')
                if warning:
                    title = compound[code[i]]['Metadata'].get('title')
                    url = compound[code[i]].get('url')
                    na_products.append({'idx' : i, 'code' : code[i], 'title' : title, 'warning' : warning, 'url' : url})
                else:
                    pass

            elif compound[code[i]]['Available_products'] == None:
                side_note = compound[code[i]].get('Side Note')
                note = side_note.get('WARNING')
                maybe_available = side_note.get('Note')
                if maybe_available:
                    # This product may be packaged on demand. Need to contact company directly for information.
                    products.append({'idx' : i,
                                     'code' : code[i],
                                     'amount' : None,
                                     'unit' : None,
                                     'price' : None,
                                     'stock_japan' : None,
                                     'aprox_delivery_time' : None,
                                     'special_remarks' : maybe_available})

                # Next are products that are not available to purchase
                elif ('discontinued' in note) or ('現在お客様の国では販売されていません' in note):
                    title = compound[code[i]]['Metadata'].get('title')
                    url = compound[code[i]].get('url')
                    na_products.append({'idx' : i, 'code' : code[i], 'title' : title, 'warning' : note, 'url' : url})

                elif ('not available for purchase' in note):
                    title = compound[code[i]]['Metadata'].get('title')
                    url = compound[code[i]].get('url')
                    na_products.append({'idx' : i, 'code' : code[i], 'title' : title, 'warning' : note, 'url' : url})

                elif ('not be available in Japan' in note):
                    title = compound[code[i]]['Metadata'].get('title')
                    url = compound[code[i]].get('url')
                    na_products.append({'idx' : i, 'code' : code[i], 'title' : title, 'warning' : note, 'url' : url})

                elif ('現在、価格および在庫状況を閲覧できません' in note) or ('not currently available' in note):
                    # Need to contact company directly for information.
                    products.append({'idx' : i,
                                     'code' : code[i],
                                     'amount' : None,
                                     'unit' : None,
                                     'price' : None,
                                     'stock_japan' : None,
                                     'aprox_delivery_time' : None,
                                     'special_remarks' : 'Price and availability are currently unavailable'})
                else:
                    print('Check compound: {}'.format(compound[code[i]]))

            else:
                for product in compound[code[i]]['Available_products']:

                    # Check if there are any special notes on the product
                    remarks = product.get('note')
                    if remarks and (remarks[-1] == '0') and (',' in remarks):
                        note = None
                    else:
                        note = remarks

                    # Get the quatity and the units of the product.
                    quantity = re.findall(r'\d+',product['ID - quatity'])
                    if quantity and len(quantity) <= 5:
                        amount = quantity[-1]
                        unit = product['ID - quatity'].split(amount)[-1].strip()
                    else:
                        amount, unit = None, None

                    # Check if there are products that need to be shipped from USA
                    delivery_time = product.get('shipment')
                    if '申し訳ございませんが、この製品のフルフィルメントと配送が遅延しています。' in delivery_time:
                        delivery_time = 'Delayed'
                    elif ('Orders outside of US' in delivery_time) or ('米国および欧州外の注文の場合' in delivery_time):
                        stock_japan = False
                        delivery_time = 'Shipped from the USA, may take several weeks.'
                    elif delivery_time:
                        stock_japan = True
                    else:
                        stock_japan = False

                    # Price might not be available and the company might need to be contacted.
                    if (product.get('price') == 'お問い合わせ') or not product.get('price'):
                        price = None
                    else:
                        price = int(product['price'].replace(',',''))

                    # Error from the scraping algorithm, same product, same units and amount mixed.
                    if (i>0) and (delivery_time[-1] == '0') and (',' in delivery_time) and not products[-1].get('price'):
                        #print('last product: ', products[-1])
                        #print('current product: ', product)
                        products[-1]['price'] = int(delivery_time.replace(',','')[1:])
                        #print('last -1 price: ', products[-1]['price'])

                    else:
                        note_warning = compound[code[i]].get('Side Note')
                        if note_warning and (('discontinued' in note_warning) or \
                                             ('現在お客様の国では販売されていません' in note_warning)):
                            pass
                        else:
                            products.append({'idx' : i,
                                            'code' : code[i],
                                            'amount' : amount,
                                            'unit' : unit,
                                            'price' : price,
                                            'stock_japan' : stock_japan,
                                            'aprox_delivery_time' : delivery_time,
                                            'special_remarks' : note})

        except ValueError:
            print(i, code[i])
            print(f'WARNING ValueError: Check compound {code[i]}, ', compound)
            pass
        except TypeError:
            print(i, code[i])
            print(f'WARNING TypeError: Check compound {code[i]}, ', compound)
            pass
        except IndexError:
            print(i, code[i])
            print(f'WARNING IndexError: Check compound {code[i]}, ', compound)
            pass
        except KeyError:
            print(i, code[i])
            print(f'WARNING KeyError: Check compound {code[i]}, ', compound)
            pass

    products_df = pd.DataFrame(products)
    na_products_df = pd.DataFrame(na_products)

    #set the type of columns to manage memory efficiently
    products_df['code'] = products_df['code'].astype('string')
    products_df['unit'] = products_df['unit'].astype('category')
    products_df['stock_japan'] = products_df['stock_japan'].astype('bool')
    return products_df.drop_duplicates(), na_products_df.drop_duplicates()
