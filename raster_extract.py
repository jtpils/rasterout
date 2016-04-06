import pandas as pd
import numpy as np
import berrl as bl
from PIL import Image
import itertools
import time
import os
import geohash


# reads the metadata txt file into memory and splits into constituent lines
def read_meta(metafile):
	with open(metafile,'rb') as f:
		f=f.read()
		f= str.split(f,'\n')
	return f 

# returns a list of the corner points of the image from teh meta data txt file
def get_corners(meta_lines):
	corner_points = [['LONG','LAT','POS']]
	for row in meta_lines:
		if 'CORNER_UL_LAT_PRODUCT' in row:
			row = str.split(row,' ')
			ul_lat = float(row[-1])
		elif 'CORNER_UL_LON_PRODUCT' in row:
			row = str.split(row,' ')
			ul_long = float(row[-1])
		elif 'CORNER_UR_LAT_PRODUCT' in row:
			row = str.split(row,' ')
			ur_lat = float(row[-1])
		elif 'CORNER_UR_LON_PRODUCT' in row:
			row = str.split(row,' ')
			ur_long = float(row[-1])
		elif 'CORNER_LL_LAT_PRODUCT' in row:
			row = str.split(row,' ')
			ll_lat = float(row[-1])
		elif 'CORNER_LL_LON_PRODUCT' in row:
			row = str.split(row,' ')
			ll_long = float(row[-1])
		elif 'CORNER_LR_LAT_PRODUCT' in row:
			row = str.split(row,' ')
			lr_lat = float(row[-1])
		elif 'CORNER_LR_LON_PRODUCT' in row:
			row = str.split(row,' ')
			lr_long = float(row[-1])

	# upper left corner 
	ul_point = [ul_long,ul_lat,'ul']
	corner_points.append(ul_point)

	# upper right corner 
	ur_point = [ur_long,ur_lat,'ur']
	corner_points.append(ur_point)

	# lower left corner 
	ll_point = [ll_long,ll_lat,'ll']
	corner_points.append(ll_point)

	# lower right corner 
	lr_point = [lr_long,lr_lat,'lr']
	corner_points.append(lr_point)


	return corner_points

# reads metadata file to memory and returns the corners
def generate_corners(dir):
	metafile = get_metafile(dir)
	meta_lines = read_meta(metafile)
	corners = get_corners(meta_lines)
	return corners

# returns a set of points that traverse the linear line between two points
def generate_points_wedges(number_of_points,point1,point2):
	# getting x points
	x1,x2 = point1[0],point2[0]
	xdelta = (float(x2) - float(x1)) / float(number_of_points)
	xcurrent = x1

	# getting y points
	y1,y2 = point1[1],point2[1]
	ydelta = (float(y2) - float(y1)) / float(number_of_points)
	ycurrent = y1

	newlist = [['LONG','LAT'],[x1,y1]]

	count = 0
	while count < number_of_points:
		count += 1
		xcurrent += xdelta
		ycurrent += ydelta
		newlist.append([xcurrent,ycurrent])

	newlist.append([x2,y2])

	return newlist

# getting lat and long for each point
def getlatlong(header,row):
	count=0
	oldrow=row
	for row in header:
		if 'lat' in str(row).lower():
			latpos = count
		elif 'long' in str(row).lower():
			longpos = count
		count+=1

	lat = float(oldrow[latpos])
	long = float(oldrow[longpos])

	return [long,lat]

# generating bottom and top points between left and right corners
# for top and bottum respectively 
def generate_horizontal_ranges(corners,widthsize):
	# corners is the output from get_corners()
	# width size is an integer for the size of the image
	
	# getting points header
	header = corners[0]

	# getting upper left point 
	upperleft  = getlatlong(header,corners[1])

	# getting upper right 
	upperright = getlatlong(header,corners[2])

	# getting lower left point 
	lowerleft = getlatlong(header,corners[-2])

	# getting lower left point 
	lowerright = getlatlong(header,corners[-1])

	# getting the upper points with the edges 
	upperpoints = generate_points_wedges(widthsize-2,upperleft,upperright)

	# getting bottom points 
	bottompoints = generate_points_wedges(widthsize-2,lowerleft,lowerright)

	return [upperpoints,bottompoints]


# generating points inbetween a set of top and bottom points 
def generate_vertical_ranges(toppoint,bottompoint,heightsize):
	points = generate_points_wedges(heightsize-2,toppoint,bottompoint)
	return points



# gets all images in a directory 
# except image b8 as its a different resolution
def get_images(dir):
	images = bl.get_filetype(dir,'TIF')
	new_images = []
	for row in images:
		#if '_B1.TIF' in str(row) or '_B2.TIF' in str(row) or '_B3.TIF' in str(row) or '_B4.TIF' in str(row):
		if '_bands_' in row:
			new_images.append(row)
	return new_images

# from a list of images returns the bands index in a list for each band possible
# this list will be added to the header value
def generate_bands(images):
	bands = []
	for row in images:
		row = str.split(row,'.')
		row = str.split(row[0],'_')
		bands.append(row[-1])
	return bands

# generator 
def gener(list):
	for row in list:
		yield row



# gets the text file that contains meta data from an output landsat directory
def get_metafile(dir):
	filename = bl.get_filetype(dir,'txt')
	filename = filename[0]
	return filename

# reads a tif file into memory and returns a dataframe of the array 
# with the appropriate index and column labels 
def generate_imageframe_single(filename):
	# reading tif image into memory and turning into numpy array
	im = Image.open(filename)
	imarray = np.array(im)
	# getting dimmensions of the image array
	dims =  imarray.shape

	# getting the size of each column and indices size
	indexsize = dims[0]
	columnsize = dims[1]

	# taking the np array to a dataframe
	imageframe = pd.DataFrame(imarray)

	# creating ranges that will be columns and indicies
	newindex = range(0,indexsize)
	newcolumns = range(0,columnsize)

	# setting the columns and indices to the lists just generated
	imageframe.columns = newcolumns
	imageframe.index = newindex

	return imageframe

# reads a tif file into memory and returns a dataframe of the array 
# with the appropriate index and column labels 
def generate_imageframes_rgb(filename):
	imageframes = []
	print filename

	if len(filename) == 1:
		filename = filename[0]



	# reading tif image into memory and turning into numpy array
	im = Image.open(filename)
	imarray = np.array(im)
	
	# splitting the array into 3 constituent np array 
	# representing colors colors blue,green,red
	imarray = np.dsplit(imarray, 3)
	for row in imarray:
		dims = row.shape
		row = row.reshape(dims[0],dims[1])


		# getting dimmensions of the image array
		dims =  row.shape

		# getting the size of each column and indices size
		indexsize = dims[0]
		columnsize = dims[1]

		# taking the np array to a dataframe
		imageframe = pd.DataFrame(row)

		# creating ranges that will be columns and indicies
		newindex = range(0,indexsize)
		newcolumns = range(0,columnsize)

		# setting the columns and indices to the lists just generated
		imageframe.columns = newcolumns
		imageframe.index = newindex

		# appending image frame to imageframe list 
		imageframes.append(imageframe)

	return imageframes

# with the appropriate index and column labels 
def generate_rgb_array(filename):
	imageframes = []

	if len(filename) == 1:
		filename = filename[0]



	# reading tif image into memory and turning into numpy array
	im = Image.open(filename)
	imarray = np.array(im)	

	return imarray

# given a point x,y in an array and a dataframe returns a intensity value
def get_point(data,x,y):
	data = data[y:y+1]
	data = data.values.tolist()
	data = data[0]
	data = data[x]
	
	return data

def get_point2(data,x,y):
	data = data[y]
	data = data[x]
	return data
# given a list of image files generates a list of dataframes all of the same size/dimmensions
def generate_all_imageframes(image_files):
	image_frames = []
	for row in image_files:
		print row
		frame = generate_imageframe(row)
		image_frames.append(frame)
	return image_frames

# given a list of image frames and a point in an x/y positon in the frame
# returns a list of intensities for each frame 
def generate_intensities(imageframes,x,y):
	intensities = []
	for row in imageframes:
		row = get_point(row,x,y)
		intensities.append(row)
	return intensities

# the get intensities function for a  3-dimmensionsal np array object 
def generate_intensities3(array,x,y):
	intensities = array[y,x]
	return intensities

# gets the rgb color band csv file for a pre-allocated tif band image file
def generate_band(folder):
	# getting image filenames 
	images = get_images(folder)

	# generating bands
	bands = generate_bands(images)

	# generating imageframes
	imageframe = generate_rgb_array(images)

	# generating shape ranges (previously columns and indexs)
	dims = imageframe.shape
	datacolumns = range(0,dims[1])
	dataindex = range(0,dims[0])


	# setting up newlists header
	header = ['X','Y','LONG','LAT','GEOHASH','RED','GREEN','BLUE']
	newlist = [header]

	# getting corner points
	corners = generate_corners(folder)

	# getting iterables of top and bottom points
	toppoints,bottompoints = generate_horizontal_ranges(corners,len(datacolumns))

	# creating generators for toppoints, bottompoints, and datacolumns
	genertop = gener(toppoints[1:])
	generbottom = gener(bottompoints[1:])
	genercolumns = gener(datacolumns)

	start = time.time()
	toppoints = []
	bottompoints =[]
	datacolumns = []

	indouter = 0 
	while indouter == 0:
		try:
			# getting row position and setting up generator
			x = next(genercolumns)
			toppoint = next(genertop)
			bottompoint = next(generbottom)

			# setting up generator for index
			generx = gener(dataindex)

			# getting points horizontal and setting up generator
			pointsindex = generate_vertical_ranges(toppoint,bottompoint,len(dataindex))

			# setting up generator to iterate through the y values
			generpoints = gener(pointsindex[1:])

			# getting the size of points index
			indexsize = len(pointsindex[1:])
			pointsindex = []

			ind = 0
			while ind == 0:
				try:
					y = next(generx)
					point = next(generpoints)
					hash = geohash.encode(float(point[1]), float(point[0]),7)
					values = generate_intensities3(imageframe,x,y).tolist()
					newlist.append([x,y,point[0],point[1],hash]+values)
				except StopIteration:
					ind = 1
					print '[%s/%s]' % (x,indexsize)
		except StopIteration:
			indouter = 1


	newlist=bl.list2df(newlist)
	newlist = newlist[(newlist.BLUE > 0)|(newlist.RED > 0)|(newlist.GREEN > 0)]
	newlist.to_csv('point_band_tiff.csv',index = False)

# gets mins from a list
def get_min(list):
	minimum = 1000
	for row in list:
		if float(row) < minimum:
			minimum = float(row)
	return minimum

# gets max from a list
def get_max(list):
	maximum = -1000
	for row in list:
		if float(row) > maximum:
			maximum = float(row)
	return maximum

# gets extrema from a geohash and table of geohashs
def get_extrema_geohash(geohash,data):
	data = data[data.GEOHASH == geohash]
	header = data.columns.values.tolist()
	data = data.values.tolist()

	lats = []
	longs = []
	for a,b in itertools.izip(header,data[0]):
		if 'LAT' in a:
			lats.append(b)
		elif 'LONG' in a:
			longs.append(b)


	latmin,latmax = get_min(lats),get_max(lats)
	longmin,longmax = get_min(longs),get_max(longs)

	extremadict = {'n':latmax,'s':latmin,'w':longmin,'e':longmax}

	return extremadict


# generate equivalent 
def generate_equilivalent(value,dimmension,min,max):
	top = value - min 
	bottom = max - min
	ratio = float(top)/float(bottom)
	value = ratio * dimmension
	value = round(value,0)
	return int(value)


# generate equivalent for y values where they start at 0 and go down
def generate_equilivalent2(value,dimmension,min,max,intialmax):
	top = value - min
	bottom = max - min
	ratio = float(top)/float(bottom)
	value = (1-ratio) * dimmension
	value = round(value,0)
	return int(value)

# given a rbg filename and an extrema returns corresponding rows and points for 
# for values that lie within extrema for the image in question
def generate_point_range(dims,folder,extremadict):
	# getting corners 
	corners = generate_corners(folder)
	initialmax = corners[1]
	initialmax = initialmax[1]

	# getting min longs  and translating column values
	corners = pd.DataFrame(corners[1:],columns = corners[0])
	longmin,longmax = corners['LONG'].min(),corners['LONG'].max()
	dimmensionx = dims[1]

	# getting range that extrema dict falls within
	x1 = generate_equilivalent(extremadict['w'],dimmensionx,longmin,longmax)
	x2 = generate_equilivalent(extremadict['e'],dimmensionx,longmin,longmax)

	# getting min latsand translating to index values
	latmin,latmax = corners['LAT'].min(),corners['LAT'].max()
	dimmensiony = dims[0]

	y1 = generate_equilivalent2(extremadict['n'],dimmensiony,latmin,latmax,initialmax)
	y2 = generate_equilivalent2(extremadict['s'],dimmensiony,latmin,latmax,initialmax)

	return [x1,x2,y1,y2]

# gets the rgb color band csv file for a pre-allocated tif band image file
def generate_band_extrema(folder,extrema):
	# getting image filenames 
	images = get_images(folder)

	# generating bands
	bands = generate_bands(images)

	# generating imageframes
	imageframe = generate_rgb_array(images)

	# generating shape ranges (previously columns and indexs)
	dims = imageframe.shape
	datacolumns = range(0,dims[1])
	dataindex = range(0,dims[0])


	# generating indices for the extrema given
	x1,x2,y1,y2 = generate_point_range(dims,folder,extrema)
	print x1,x2,y1,y2

	# setting up newlists header
	header = ['X','Y','LONG','LAT','GEOHASH','RED','GREEN','BLUE']
	newlist = [header]

	# getting corner points
	corners = generate_corners(folder)

	# getting iterables of top and bottom points
	toppoints,bottompoints = generate_horizontal_ranges(corners,len(datacolumns))

	# creating generators for toppoints, bottompoints, and datacolumns
	genertop = gener(toppoints[1+x1:1+x2])
	generbottom = gener(bottompoints[1+x1:1+x2])
	genercolumns = gener(datacolumns[x1:x2])

	start = time.time()
	toppoints = []
	bottompoints =[]
	datacolumns = []

	indouter = 0 
	while indouter == 0:
		try:
			# getting row position and setting up generator
			x = next(genercolumns)
			toppoint = next(genertop)
			bottompoint = next(generbottom)

			# setting up generator for index
			generx = gener(dataindex)

			# getting points horizontal and setting up generator
			pointsindex = generate_vertical_ranges(toppoint,bottompoint,len(dataindex))

			# setting up generator to iterate through the y values
			generpoints = gener(pointsindex[1+y1:1+y2])

			# getting the size of points index
			indexsize = len(pointsindex[1+y1:1+y2])
			pointsindex = []

			ind = 0
			while ind == 0:
				try:
					y = next(generx)
					point = next(generpoints)
					hash = geohash.encode(float(point[1]), float(point[0]),7)
					values = generate_intensities3(imageframe,x,y).tolist()
					newlist.append([x,y,point[0],point[1],hash]+values)
				except StopIteration:
					ind = 1
					print '[%s/%s]' % (x,indexsize)
		except StopIteration:
			indouter = 1


	newlist=bl.list2df(newlist)
	print newlist
	newlist = newlist[(newlist.BLUE > 0)|(newlist.RED > 0)|(newlist.GREEN > 0)]
	newlist.to_csv('point_band_tiff.csv',index = False)
	return newlist