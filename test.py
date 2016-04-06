import pandas as pd
import numpy as np
import berrl as bl
from raster_extract import *
key='pk.eyJ1IjoibXVycGh5MjE0IiwiYSI6ImNpam5kb3puZzAwZ2l0aG01ZW1uMTRjbnoifQ.5Znb4MArp7v3Wwrn6WFE6A'


# reading square csv file into memory 
squares = pd.read_csv('squares5.csv')


# getting extrema dictionary to send into the generate band function
extremadict = get_extrema_geohash('dnwke',squares)

# getting the geohash dataframe so it can be written to a block
block = squares[squares.GEOHASH == 'dnwke']

# generating points
#extremadict = {'s': 37.353515625, 'e': -81.03515625, 'w': -81.0791015625, 'n': 37.3974609375}

dims = (15661, 15381, 3)

data = generate_band_extrema('meta',extremadict)
bl.make_points(data,list=True,filename='corners.geojson')
bl.make_blocks(block,list=True,filename='block.geojson')
bl.loadparsehtml(bl.collect(),key)