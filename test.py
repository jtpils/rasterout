import berrl as bl
import pandas as pd
from pandas_sparse import *
import numpy as np

# making a dataframe of all pixels in the image
mapdf = make_dfmap('b.TIF')
print mapdf

# making all the pixels into a geojson file
bl.make_blocks(mapdf[:30000],list=True,filename='blocks.geojson',bounds=True)

# making a file dictionary object used to parse the html/js
file_dictionary = bl.make_file_dict('blocks.geojson',mapdf,['GEOHASH'],slider_fields=['RED','GREEN','BLUE'],zooms=[10,20])

# loading map based on inputs 
bl.loadparsehtml(bl.collect(),True,file_dictionary=file_dictionary,sidebar=True,bounds=True,colorkey='COLORKEY')



