import pandas as pd
import numpy as np
import berrl as bl
import raster_extract



# making image dataframe
data = raster_extract.make_image('b.tif')

# taking the red,green,blue fields and turing them into a color hash
data = bl.make_colorkey_table(data)

# making the geometries out of the point data
data = raster_extract.make_squares_table(data)

bl.make_blocks(data,list=True,filename='blocks.geojson')

bl.loadparsehtml(['blocks.geojson'],True,colorkey=True)

