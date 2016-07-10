from raster_extract import *
import numpy as np
import pandas as pd
import berrl as bl
import pyproj
from math import radians, degrees
import math
import rasterio
import warnings
warnings.simplefilter(action = "ignore", category = Warning)
import time

# gets the distance between two points
def distance2points(point1,point2):
	# getting point1 lat andlong
	long1,lat1 = [point1[0],point1[1]]

	# getting point1 lat andlong
	long2,lat2 = [point2[0],point2[1]]
	
	# calculating delta theta
	theta = lat2 - lat1 
	theta = radians(theta)

	# calculating delta sigma
	sig = long2 - long1
	sig = radians(sig)

	# converting to radians
	lat1 = radians(lat1)
	long1 =radians(long1)
	lat2 = radians(lat2)
	long2 = radians(long2)

	# radius in meters
	r = 6371000.0 # meters

	a = (math.sin(theta/2) ** 2) + (math.cos(lat1) * math.cos(lat2) * math.sin(sig/2)**2)

	c = 2 * math.atan2(a ** .5, (1 - a)**.5)

	distance = r * c
	return distance


# given f ratio (the array positon betweeen the top and bottom point) divide by the total len() of the column
# returns the x and y latitude and longitude
def intermediate_point(f,sigma,point1,point2):
	# getting point1 lat andlong
	long1,lat1 = [point1[0],point1[1]]

	# getting point1 lat andlong
	long2,lat2 = [point2[0],point2[1]]

	# converting to radians
	lat1 = radians(lat1)
	long1 =radians(long1)
	lat2 = radians(lat2)
	long2 = radians(long2)

	# finding a  
	a = math.sin((1 - f) * sigma) / math.sin(sigma)

	# finding b 
	b = math.sin(f * sigma) / math.sin(sigma)


	x = (a * math.cos(lat1) * math.cos(long1)) + (b * math.cos(lat2) * math.cos(long2))
	y = (a * math.cos(lat1) * math.sin(long1)) + (b * math.cos(lat2) * math.sin(long2))
	z = (a * math.sin(lat1)) + (b * math.sin(lat2))

	lat = math.atan2(z,(x**2 + y**2)**.5)
	long = math.atan2(y,x)

	x,y = degrees(long),degrees(lat)

	return [x,y]


# generating a set of points from the same inputs as generate_point_wedges
# attempting to project differently and still function the same
def generate_points_wedges2(number_of_points,point1,point2,**kwargs):
	dimmension = False
	bb = False
	for key,value in kwargs.iteritems():
		if key == 'dimmension':
			dimmension = value
		if key == 'bb':
			bb = value

	# getting distance
	distance = distance2points(point1,point2)

	# radius of the earthin meters
	r = 6371000.0 # meters

	# finding sigma
	sigma = float(distance) / r

	# setting up newlist and count
	count = 0
	newlist = [['LONG','LAT']]

	if not dimmension == False and not bb == False:
		if dimmension == 'x':
			newlist.append(point1[0])
		elif dimmension == 'y':
			newlist.append(point1[1])

	
	# iterating through each point calculating the percentage of pixils or nparray rows completed
	# as the percentage of circle traversed f
	while not len(newlist) == number_of_points:
		f = float(count) / float(number_of_points)
		#print f,sigma
		#raw_input()
		point = intermediate_point(f,sigma,point1,point2)
		count += 1

		if not dimmension == False:
			if dimmension == 'y':
				newlist.append(point[1])
			elif dimmension == 'x':
				newlist.append(point[0])
			else:
				newlist.append(point)

	if not dimmension == False:
		if bb == False:
			if dimmension == 'x':
				newlist.append(point2[0])
			elif dimmension == 'y':
				newlist.append(point2[1])
		return newlist[1:]
	else:
		return newlist


# makes the master list of indicies that will be used to build the entire table
def make_indicies(x,y):
	xs = range(x)
	ys = range(y)
	newlist = []
	for row in xs:
		row = list(zip([row]*len(ys),ys))
		row = pd.DataFrame(row,columns=['X','Y'])
		row['string'] = row['X'].astype(str) + ',' + row['Y'].astype(str)
		row = row['string'].values.tolist()
		newlist.append(row)

	newlist = pd.DataFrame(newlist)
	newlist = newlist.T
	return newlist

def make_points(x,y,extrema):
	newlist = []
	extremaoppositex = [extrema['n'],extrema['s']]
	extremaoppositey = [extrema['w'],extrema['e']]

	xs1 = generate_points_wedges2(x,[extrema['w'],40.0],[extrema['e'],40.0],dimmension='x')
	ys1 = generate_points_wedges2(y,[-80,extrema['n']],[-80,extrema['s']],dimmension='y')

	for row in xs1:
		row = list(zip([row]*len(ys1),ys1))
		row = pd.DataFrame(row,columns=['X','Y'])
		row['string'] = row['Y'].astype(str) + ',' + row['X'].astype(str)
		row = row['string'].values.tolist()

		newlist.append(row)
	newlist = pd.DataFrame(newlist)
	newlist = newlist.T

	return newlist

# this function is quite literally where the entire mapped dataframe is created
def map_src(arg):
	global datasrc
	global pointsdf
	argstring = arg
	# getting out indices to be used in other inputs for getting rgb/point
	arg = str.split(arg,",")
	x,y = int(arg[0]),int(arg[1])

	# 
	rgb = datasrc[:,y,x]
	rgb = str(rgb[0]) + ',' + str(rgb[1]) + ',' + str(rgb[2])

	# getting the point value of the center point
	point = pointsdf[x][y]

	stringoverwrite = argstring + ',' + rgb + ',' + point 

	return stringoverwrite


# given a source frame returns the delta of talues based on the pixel spacing
def get_range_delta(pointdf):
	# getting the 2 points to compare
	point1 = pointdf.iloc[0][0]
	point2 = pointdf.iloc[1][1]

	# getting the actual point values instead of the stings
	point1 = [float(str.split(point1,',')[-1]),float(str.split(point1,',')[-2])]
	point2 = [float(str.split(point2,',')[-1]),float(str.split(point2,',')[-2])]

	distance_vert = abs(point1[1] - point2[1]) / 2.0
	distance_horz  = abs(point1[0] - point2[0]) / 2.0

	return distance_vert,distance_horz


# this function is mapped against the string df and yields a 
# series containing every block value found in the image
def split(each):
	each = str.split(each,',')
	return [int(each[0]),int(each[1]),float(each[2]),float(each[3]),int(each[4]),int(each[5]),int(each[6])]


# returns a dfmap of every pixel in in an a geospatial raster image
# this can be sent directly into make_blocks and already contains a colorfield
def make_dfmap(filename,**kwargs):
	start = time.time()
	dims = False
	extrema = False
	for key,value in kwargs.iteritems():
		if key == 'extrema':
			extrema = value
		if key == 'dims':
			dims = value

	# reading the rasterio datasrc into memory
	datasrc = rasterio.open(filename)

	# assigning dims if not set by kwarg
	if dims == False:	
		dims = datasrc.shape


	# asigning extrema if not set by kwarg
	if extrema == False:
		if 'epsg:4326'==str(datasrc.meta['crs']['init']):
			# getting extrema
			extrema = get_extrema(datasrc.bounds)
		else:
			extrema = get_extrema(datasrc.bounds,transform=True)

	#extrema = get_extrema(datasrc.bounds)
	global r,g,b
	global pointsdf
	# reading red,green,blue values into memory
	r, g, b =  datasrc.read()

	# reshaping each r,g,b value to a flat list
	numberx = dims[1]
	numbery = dims[0]
	#r,g,b = np.array(r),np.array(g),np.array(b)
	#r,g,b = r.reshape(r.shape[0]*r.shape[1]),g.reshape(g.shape[0]*g.shape[1]),b.reshape(b.shape[0]*b.shape[1])
	
	# the points df will be indexed to get points into the the index df which will have the apply_map shit done
	pointsdf = make_points(numberx,numbery,extrema)
	distance_vert,distance_horz = get_range_delta(pointsdf) 
	#print time.time() - start,'pointsdf'

	#twogroups = '(\,+)(?P<digit>[0-10000])(?P<digit>[0-10000])'
	#print re.split('\,+')
	# getting np matrix 
	inddf = make_indicies(numberx,numbery)
	inddf = inddf.unstack(level=0)
	#print time.time() - start,'inddf'
	#print inddf.shape
	#print inddf + ',' + pointsdf.unstack(level=0)

	# making string/ind dataframe that will be applied to rgb values
	df = pd.DataFrame({'ind':inddf,'string':inddf + ',' + pointsdf.unstack(level=0)})
	df = df.reset_index()

	# expanding out values
	expanded = df['ind'].str.split(',',expand=True)
	#print time.time() - start,'expanded'

	#print df.index.str.extract(twogroups, expand=True)
	df[['x','y']] = expanded.astype(int)
	#print time.time() - start,'astype_int'
	
	# mapping rbg values about flat index
	
	df['string'] = df['string'] + ',' + r[df['y'],df['x']].astype(str) + ',' + g[df['y'],df['x']].astype(str) + ',' + b[df['y'],df['x']].astype(str)
	#print time.time() - start,'rgb'
	
	# splitting the string field into list subsets 
	# replace this later
	newdf = df['string'].str.split(',',expand=True)
	newdf = newdf.astype(float) 
	df = pd.DataFrame(newdf.values.tolist(),columns=['X','Y','LAT','LONG','RED','GREEN','BLUE'])
	#print time.time() - start,'final_split'


	# syntax se,sw,ne,nw
	df['LAT1'] = df['LAT'] - distance_vert
	df['LONG1'] = df['LONG'] + distance_horz
	df['LAT2'] = df['LAT'] - distance_vert
	df['LONG2'] = df['LONG'] - distance_horz
	df['LAT3'] = df['LAT'] + distance_vert
	df['LONG3'] = df['LONG'] + distance_horz
	df['LAT4'] = df['LAT'] + distance_vert
	df['LONG4'] = df['LONG'] - distance_horz
	#print time.time() - start,'points'

	# geohashing the dataframe against the lat and long fields
	df = bl.map_table(df,8)

	# making colorkey table
	df = bl.make_colorkey_table(df)

	# making colors integers
	df[['RED','GREEN','BLUE']] = df[['RED','GREEN','BLUE']].astype(int)


	headers = []
	for row in df.columns.values.tolist():
		if not 'LAT' == str(row) and not 'LONG' == str(row):
			headers.append(row)

	return df[headers]




def make(inddf,pointsdf):
	pts_splits = np.array_split(pointsdf,16)
	ind_splits = np.array_split(inddf,16)
	return zip(ind_splits,pts_splits)

def split_dfmap(filename):
	datasrc = rasterio.open(filename)

	dims = datasrc.shape

	if 'espg:4326' in str(datasrc.meta['crs']['init']):
		# getting extrema
		extrema = get_extrema(datasrc.bounds)
	else:
		extrema = get_extrema(datasrc.bounds,transform=True)
	#extrema = get_extrema(datasrc.bounds)

	datasrc = datasrc.read()

	numberx = dims[1]
	numbery = dims[0]

	# the points df will be indexed to get points into the the index df which will have the apply_map shit done
	pointsdf = make_points(numberx,numbery,extrema)
	inddf = make_indicies(numberx,numbery)

	global datasrc

	return make(inddf,pointsdf)

# previou
def mapped(args):
	inddf, pointsdf,filename = args
	global pointsdf
	dfmap = inddf.applymap(map_src)
 	#expand_dfmap(dfmap).to_csv(str(random.randint(0,1000)) + '.csv',index=False)
 	dfmap.to_csv(filename,index=False)
 	return []


# given a x or y min or max and a pixel size returns a definitive range <= >= list of the values 
# that will sliced in each dimmension iteration
def make_ranges(min,max,pixelsize):
	rangelist = []
	current = min
	oldcurrent = min
	while current <= max:
		current += pixelsize
		rangelist.append([oldcurrent,current-1])
		oldcurrent = current

	rangelist = rangelist[:-1] + [[rangelist[-2][1],max]]

	return rangelist

# returns a slice of dataframe 
# this will be the first slice about an axis 
# hopefully using a smaller df for each individual slice will help
def slicey(miny,maxy,data):
	partialdata = data[(data.Y >= miny)&(data.Y <= maxy)]
	return partialdata

# returns a slice of dataframe 
# this will be the first slice about an axis 
# hopefully using a smaller df for each individual slice will help
def slicex(minx,maxx,partialdf):
	downsampled_df = partialdf[(partialdf.X >= minx)&(partialdf.X <= maxx)]
	return downsampled_df

# previously written function for downsampling image pixels give fields
def downsample(data,pixelsize):
	miny,maxy = data.Y.min(),data.Y.max()
	minx,maxx = data.X.min(),data.X.max()

	yranges = make_ranges(miny,maxy,pixelsize)
	xranges = make_ranges(minx,maxx,pixelsize)

	header = data.columns.values.tolist()
	# getting newheaders out of initial header
	newlist = []
	for row in header:
		if not row == 'LAT' and not row =='LONG':
			newlist.append(row)

	data = data[newlist]


	downsampled_list = [header]

	x = 0 
	y = 0
	for row in yranges:
		oldrow = row
		# slicing by the y component to yield a smaller df
		partialdf = slicey(row[0],row[1],data)
		y += 1
		x = 0
		for row in xranges:
			downsampled_df = slicex(row[0],row[1],partialdf)

			extrema = {'s':downsampled_df['LAT1'].min(),'w':downsampled_df['LONG1'].min(),'n':downsampled_df['LAT3'].max(),'e':downsampled_df['LONG4'].max()}
			points = [extrema['s'],extrema['w'],extrema['s'],extrema['e'],extrema['n'],extrema['w'],extrema['n'],extrema['e']]
			geohash = downsampled_df['GEOHASH'][:1].values.tolist()[0][:-2]
  
			red,green,blue = int(downsampled_df['RED'].mean()),int(downsampled_df['GREEN'].mean()),int(downsampled_df['BLUE'].mean())
			colorkey = bl.hexify(red,green,blue)

			values = [x,y,red,green,blue,colorkey]
			newrow = [geohash] + points + values
			downsampled_list.append(newrow)



	downsampled_list = bl.list2df(downsampled_list)
	#downsampled_list.columns = newlist

	return downsampled_list

# map a an iteration of a down sample process meant to be parrellized 
# in spark
def map_downsample(args):
	data,yranges,y,pixelsize,filename = args
	header = data.columns.values.tolist()
	minx,maxx = data.X.min(),data.X.max()
	xranges = make_ranges(minx,maxx,pixelsize)
	header = ['GEOHASH1'] + header[1:]
	downsampled_list = [header]

	for row in yranges:
		oldrow = row
		# slicing by the y component to yield a smaller df
		partialdf = slicey(row[0],row[1],data)
		y += 1
		x = 0
		for row in xranges:
			downsampled_df = slicex(row[0],row[1],partialdf)
			if len(downsampled_df) > 0:
				extrema = {'s':downsampled_df['LAT1'].min(),'w':downsampled_df['LONG1'].min(),'n':downsampled_df['LAT3'].max(),'e':downsampled_df['LONG4'].max()}

				points = [extrema['s'],extrema['w'],extrema['s'],extrema['e'],extrema['n'],extrema['w'],extrema['n'],extrema['e']]
				geohash = downsampled_df['GEOHASH'][:1].values.tolist()[0][:-2]
				x += 1



				red,green,blue = int(downsampled_df['RED'].mean()),int(downsampled_df['GREEN'].mean()),int(downsampled_df['BLUE'].mean())
				colorkey = bl.hexify(red,green,blue)

				values = [x,y,red,green,blue,colorkey]
				newrow = [geohash] + points + values
				downsampled_list.append(newrow)
	data = []
	downsampled_list = bl.list2df(downsampled_list)
	return downsampled_list

# creates the arguments in which spark will parrellize the mapping of 
# about the mapped function above
def make_downsample_args(data,pixelsize,splits):
	argslist = []
	miny,maxy = data.Y.min(),data.Y.max()
	minx,maxx = data.X.min(),data.X.max()

	yranges = make_ranges(miny,maxy,pixelsize)
	xranges = make_ranges(minx,maxx,pixelsize)

	current = 0
	oldcurrent = 0
	rangelist = []
	delta = len(yranges) / splits
	while not len(rangelist) == splits:
		current += delta
		rangelist.append([oldcurrent,current])
		oldcurrent = current

	rangelist = rangelist[:-1] + [[rangelist[-1][0],len(yranges)-1]]

	count = 0
	for row in rangelist:
		yrangestemp = yranges[row[0]:row[1]]
		miny,maxy = yrangestemp[0],yrangestemp[-1]
		miny,maxy = miny[0],maxy[1]
		datatemp = slicey(miny,maxy,data)
		args = [datatemp,yrangestemp,row[0],pixelsize,str(count)+'.geojson']
		argslist.append(args)
	return argslist

# spark wrapper function for downsampling api for images
def make_spark_downsample(data,pixelsize,splits,sc):
	args = make_downsample_args(data,pixelsize,splits)
	concurrent = sc.parallelize(args)
	data = concurrent.map(map_downsample).collect()
	data = pd.concat(data)
	return data


###############################################################################################################################																				

# rgb is equal 2,3,4
def red(arg):
	color = int(str.split(arg,',')[2])
	return color

# rgb is equal 2,3,4
def green(arg):
	color = int(str.split(arg,',')[3])
	return color

# rgb is equal 2,3,4
def blue(arg):
	color = int(str.split(arg,',')[4])
	return color

def make_colormaps(dfmap):
	redmap = dfmap.applymap(red)
	greenmap = dfmap.applymap(green)
	bluemap = dfmap.applymap(blue)
 
	return redmap,greenmap,bluemap

# returns averages of red,blue,green kernel ranges
def rgb_averages(redmap,green,bluemap,xindex,yindex):
	red = kernelaverage(redmap,xindex,yindex)
	green = kernelaverage(greenmap,xindex,yindex)
	blue = kernelaverage(bluemap,xindex,yindex)
	return int(red),int(green),int(blue)

#3,6 3,6
def kernelaverage(colormap,xindex,yindex):
	kernel = colormap[range(xindex[0],xindex[1])][yindex[0]:yindex[1]]
	kernel = kernel.stack().mean()
	return kernel

def make_point_middle(xindex,yindex,pointmap):
	kernel = colormap[range(xindex[0],xindex[1])][yindex[0]:yindex[1]]
	kernel = kernel.stack().mean()
	return kernel

# take note of rare syntax this for pipegeoson integration
def retrieve_long_lat(string):
	return [float(str.split(string,',')[0]),float(str.split(string,',')[1])]




