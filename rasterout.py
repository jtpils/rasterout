import rasterio
import numpy as np
import pandas as pd
import pyproj

# get bounds
def get_bounds(data,**kwargs):
	transform = False
	for key,value in kwargs.iteritems():
		if key == 'transform':
			transform = value
	north = data.top
	south = data.bottom
	west = data.left
	east = data.right

	extrema = {'n':north,'s':south,'w':west,'e':east}

	# newcorners
	ur_point = [extrema['e'],extrema['n']]
	ll_point = [extrema['w'],extrema['s']]


	if transform == True:
		p1 = pyproj.Proj(init='epsg:'+str(3857))
		p2 = pyproj.Proj(init='epsg:'+str(4326))

		ur_point = pyproj.transform(p1,p2,ur_point[0],ur_point[1])
		ll_point = pyproj.transform(p1,p2,ll_point[0],ll_point[1])

	extrema = {'n':ur_point[1],'s':ll_point[1],'w':ll_point[0],'e':ur_point[0]}
	return extrema

# make points and indicies functions
def make_points_grid(extrema,x,y):
	lats = np.linspace(extrema['s'],extrema['n'],y)
	longs = np.linspace(extrema['w'],extrema['e'],x)
	xs = range(x)
	ys = range(y)
	deltay = (extrema['n'] - extrema['s']) / float(y)
	deltax = (extrema['e'] - extrema['w']) / float(x)
	deltax,deltay = deltax / 2.,deltay / 2.


	sizelongs = len(longs)
	points  = []
	indicies = []
	for lat,y in zip(lats.tolist(),ys):
		indicies += zip(xs,[y]*sizelongs)
		points += zip(longs,[lat] * sizelongs)
	
	data = pd.DataFrame(points,columns=['LONG','LAT'])
	data[['X','Y']] = pd.DataFrame(indicies,columns=['X','Y'])[['X','Y']]

	return data,deltax,deltay

# make points and indicies functions
def make_points_grid_bounds(extrema,xranges,yranges,xtotal,ytotal):

	deltay = (extrema['n'] - extrema['s']) / float(ytotal)
	deltax = (extrema['e'] - extrema['w']) / float(xtotal)
	deltax,deltay = deltax / 2.,deltay / 2.

	
	print extrema,xranges
	# getting the slopes to calculate the real extrema
	s,n = np.interp([(ytotal-yranges[0]),(ytotal- yranges[-1])],[0,ytotal],[extrema['s'],extrema['n']])
	w,e = np.interp([xranges[0],xranges[-1]],[0.,xtotal],[extrema['w'],extrema['e']])
	
	extrema =  {'n':n,'s':s,'e':e,'w':w}

	lats = np.linspace(extrema['s'],extrema['n'],len(yranges))
	longs = np.linspace(extrema['w'],extrema['e'],len(xranges))


	sizelongs = len(longs)
	points  = []
	indicies = []
	for lat,y in zip(lats.tolist(),yranges):
		indicies += zip(xranges,[y]*sizelongs)
		points += zip(longs,[lat] * sizelongs)
	
	data = pd.DataFrame(points,columns=['LONG','LAT'])
	data[['X','Y']] = pd.DataFrame(indicies,columns=['X','Y'])[['X','Y']]

	return data,deltax,deltay



# function thing
def make_raster_df(filename,bounds=False,downscale_color=False):
	datasrc = rasterio.open(filename)

	if 'epsg:4326'==str(datasrc.meta['crs']['init']):
		# getting extrema
		extrema = get_bounds(datasrc.bounds)
	else:
		extrema = get_bounds(datasrc.bounds,transform=True)
	

	dims = datasrc.shape
	# reading rgb
	r, g, b =  datasrc.read()
	x = dims[1]
	y = dims[0]
	datasrc.close()
	
	# making points grid
	if bounds == False:
		data,deltax,deltay = make_points_grid(extrema,x,y)
	else:
		minx,maxx,miny,maxy = bounds
		xranges = range(minx,maxx)
		yranges = range(miny,maxy)
		data,deltax,deltay = make_points_grid_bounds(extrema,xranges,yranges,x,y)

	# adding red blue and green colors
	data[['RED','BLUE','GREEN']] = pd.DataFrame({"RED":r[data['Y'],data['X']],"BLUE":b[data['Y'],data['X']],"GREEN":g[data['Y'],data['X']]})[['RED','BLUE','GREEN']]

	# adding north south east and west fields
	data['NORTH'] = data['LAT'] + deltay
	data['SOUTH'] = data['LAT'] - deltay
	data['EAST'] = data['LONG'] + deltay
	data['WEST'] = data['LONG'] - deltay

	if downscale_color == True:
		data[['RED','GREEN','BLUE']] = data[['RED','GREEN','BLUE']] / 20
		data[['RED','GREEN','BLUE']] = data[['RED','GREEN','BLUE']].round(0) * 20

	# adding colorkey field
	newlist = []
	for r,g,b in data[['RED','BLUE','GREEN']].astype(int).values:
		if r > 255:
			r = 255
		if g > 255:
			g = 255
		if b > 255:
			b = 255
		r,g,b = str(hex(r))[2:],str(hex(g))[2:],str(hex(b))[2:]
		if len(r) == 1:
			r = '0' + r
		if len(g) == 1:
			g = '0' + g
		if len(b) == 1:
			b = '0' + b
		if len(r) == 0:
			r = '00'
		if len(g) == 0:
			g = '00'
		if len(b) == 0:
			b = '00'
		if str(r[-1]) == "L":
			r = r[:-1]
		if str(g[-1]) == "L":
			g = g[:-1]
		if str(b[-1]) == "L":
			b = b[:-1]
		newlist.append('#' + r + g + b)
	data['COLORKEY'] = newlist
	return data

# downsample the raster out dataset
def downsample(data,size):
	# setting up the xs and ys
	data['XL'] = data['X'] / size
	data['YL'] = data['Y'] / size
	data[['XL','YL']] = data[['XL','YL']].round(0)


	# grouped the new xs and ys
	grouped = data.groupby(['XL','YL'])
	
	# reevaluating columns
	mask = grouped.first()
	mask[['LONG','LAT','RED','BLUE','GREEN']] = grouped[['LONG','LAT','RED','BLUE','GREEN']].mean()
	mask[['NORTH','EAST']] = grouped[['NORTH','EAST']].max()
	mask[['SOUTH','WEST']] = grouped[['SOUTH','WEST']].min()
	

	return mask.reset_index()
'''
import mapkit as mk
mk.cln()
data = make_raster_df('b.tif',downscale_color=True)
mk.cln()
mk.make_map([data,'blocks'])
mk.a()
'''