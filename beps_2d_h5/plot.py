#ffmpeg -framerate 2 -i plots/charge%d.png -c:v libx264 -vf 'scale=-1:1080,format=yuv420p' out9.mp4
#ffmpeg -framerate 10 -pattern_type glob -i 'plots/*.png' -c:v libx264 -vf 'scale=-1:1080,format=yuv420p' out_out2.mp4
import matplotlib
matplotlib.use('Agg')
import h5py 
import matplotlib.ticker as ticker	
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
from joblib import Parallel, delayed
import numpy as np
import multiprocessing
import os
import glob
import shutil
import sys
import subprocess

def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'$ 10^{' +str(b)+'}$'

dpi = 300
ind = 0
def main():
	args= sys.argv
	print 'num cores: ' +str(multiprocessing.cpu_count())
	if len(args) > 1:
		print args
		data_params = read_input_deck(args[1])
		global_data = data_params[0]
		num_plots = global_data[0]

		folder = []
		startn = [global_data[1]]
		ndump = [global_data[2]]
		lastn = [global_data[3]]
		other = data_params[1]
		titles = [] 
		lineout = []
		log = []
		max_min_overwrite = []
		thresh_value = np.zeros((num_plots))
		max_over = []
		min_over = []
		threshold_over = []
		xmin = np.zeros((num_plots))
		xmax = np.zeros((num_plots))
		file_out = []
		colormap = []
		numx = int((lastn[0]-startn[0])/ndump[0])
		for j in xrange(num_plots):
			folder.append([0])
			file_out.append([])
			folders = (other[j][0])
			first_folder = folders[0]
			file_out[j].append(min(glob.iglob(first_folder+'/*.h5'),key=os.path.getctime)[:-9])
			titles.append(other[j][1])
			if(other[j][2]):
				log.append('y')
			else:
				log.append('n')

			if(other[j][3]):
				lineout.append('y')
			else:
				lineout.append('n')

			if(other[j][4]):
				max_min_overwrite.append('y')
			else:
				max_min_overwrite.append('n')

			min_over.append((other[j][5])[0])
			max_over.append((other[j][5])[1])
			threshold_over.append((other[j][5])[2])
			colormap.append(other[j][6])
		index_info = [startn,ndump,lastn,folder]
	else:
		rootdir = '.'
		plt.close("all")
		dirs = []
		file_list = []
		for subdirs in os.listdir(rootdir):
			dirs.append(subdirs)
		print_dirs(dirs)

		data_list = raw_input('Select Simulations to Plot (i.e. 0,1,2) : ')
		data_list = data_list.split(',')
		data_list = [int(a) for a in data_list]
		frac = []
		file_out2 = []
		for j in xrange(len(data_list)):
			file_out2.append(dirs[data_list[j]])
		print file_out2

		frac = []
		file_out = []
		data_list = []
		last_ind = []
		first_ind = []
		file_list = []
		dir_list = []
		directs = ''
		folder_point = []
		for j_ind in xrange(len(file_out2)):
			rootdir = file_out2[j_ind]
			if (j_ind == len(file_out2)-1):
				directs = directs + rootdir
			else:
				directs = directs + rootdir + ', '
			for root, subdirs, files in os.walk(rootdir):
				for j in files:
					if('0.h5' in j and root not in dir_list):
						dir_list.append(root)
						file_list.append(j)
						folder_point.append(j_ind)
						

		print 'Listing Data in '+ directs + ':'
		print_dirs(dir_list)
		num_plots = int(raw_input('Number of Subplots (1-): '))
		log = []
		lineout= []
		titles = []
		thresh_value = np.zeros((num_plots))
		xmin = np.zeros((num_plots))
		xmax = np.zeros((num_plots))
		folder = []
		max_min_overwrite = []
		max_over = []
		min_over = []
		threshold_over = []
		colormap = []
		for j_ind in xrange(num_plots):
			file_out.append([])
			frac.append([])
			data_list.append([])
			last_ind.append([])
			first_ind.append([])
			folder.append([])
			f = -1
			f2 = -1
			last = -1
			first = -1
			data_list[j_ind] = raw_input('Select Data for Plot #' + str(j_ind+1)+' (i.e. 0,1,2): ')
			data_list[j_ind] = data_list[j_ind].split(',')
			data_list[j_ind] = [int(a) for a in data_list[j_ind]]
			
			
			for j in xrange(len(data_list[j_ind])):
				new2 =glob.iglob((dir_list)[(data_list[j_ind])[j]]+'/*.h5')
				
				for v in new2:
					index = v.index('.')
					f_curr =int(v[index-6:index])
					#if(v[:index-6] not in file_out):
					#file_out.append(v[:index-6])
					#print file_out
					last = f_curr
					if(f == -1 ):
						first = f_curr
						f = f_curr
					elif(f2 == -1):
						f2 = 1
						f = f_curr - f
				first_ind[j_ind].append(first)
				last_ind[j_ind].append(f_curr)
				file_out[j_ind].append(v[:index-6])
				frac[j_ind].append(f)
				folder[j_ind].append(folder_point[(data_list[j_ind])[j]])
				f = -1
				f2 = -1
			print folder
			thresh_value[j_ind]= -np.inf
			titles.append(raw_input('Title for Plot #' + str(j_ind+1) +' :' ))
			ln_out = raw_input('Lineout on axis? (y/n):')
			lineout.append(ln_out)
			if(ln_out != 'y'):
				log.append(raw_input('Do you want to plot in log space? (y/n):'))
			else:
				log.append('n')

			if(log[j_ind] == 'y'):
				overwrite = raw_input('Overwrite min max? (y/n):')
				max_min_overwrite.append(overwrite)
				if(overwrite == 'y'):
					info_list = raw_input('Specify Axis info (min, max, linear threshold):')
					info_list = info_list.split(',')
					info_list = [float(a) for a in info_list]
					min_over.append(info_list[0])
					max_over.append(info_list[1])
					threshold_over.append(info_list[2])
				else:
					min_over.append(-np.inf)
					max_over.append(np.inf)
					threshold_over.append(0)
				colormap.append(raw_input('Specify a Colormap (def,bwr,..):'))
			else:
				overwrite = raw_input('Overwrite min max? (y/n):')
				max_min_overwrite.append(overwrite)
				if(overwrite == 'y'):
					info_list = raw_input('Specify Axis info (min,max):')
					info_list = info_list.split(',')
					info_list = [float(a) for a in info_list]
					min_over.append(info_list[0])
					max_over.append(info_list[1])
				else:
					min_over.append(-np.inf)
					max_over.append(np.inf)
				if(ln_out == 'y'):
					colormap.append('def')
				else:
					colormap.append(raw_input('Specify a Colormap (def,bwr,..):'))
				threshold_over.append(0)
				
			#data_list3 = raw_input('Specify Data bounds (xmin,xmax) : ')
				#data_list3 = data_list3.split(',')
				#data_list3 =  [float(a) for a in data_list3]
		ndump = np.zeros(len(file_out2))
		lastn = np.zeros(len(file_out2))
		firstn = np.zeros(len(file_out2))
		print frac
		for xj in xrange(num_plots):
			for jj in xrange(len(data_list[xj])):
				if((frac[xj])[jj] > ndump[folder[xj][jj]]):
					ndump[folder[xj][jj]] = (frac[xj])[jj]
				if((last_ind[xj])[jj] > lastn[folder[xj][jj]]):
					lastn[folder[xj][jj]] = (last_ind[xj])[jj]
				if((first_ind[xj])[jj] > firstn[folder[xj][jj]]):
					firstn[folder[xj][jj]] = (first_ind[xj])[jj]
				
		
		startn = np.zeros(len(file_out2))
		numx = np.inf
		for xj in xrange(len(file_out2)):
			print 'In ' + file_out2[xj] + ', ndump starts at ' + str(int(firstn[xj])) + ' and ends on ' + str(int(lastn[xj])) + ' in intervals of '+ str(int(ndump[xj])) 
			startn[xj] = int(raw_input('Select n start:'))
			ndump[xj] = int(raw_input('Select ndump (multiples of '+str(int(ndump[xj]))+ '):'))
			num_temp = int((lastn[xj]-startn[xj])/ndump[xj])
			if(numx == -1):
				numx = num_temp
			else:
				if(numx > num_temp):
					numx = num_temp
		index_info = [startn,ndump,lastn,folder]
	if not os.path.exists('plots'):
		os.makedirs('plots')
	else:
		shutil.rmtree('plots')
		os.makedirs('plots')		

	if(num_plots > 3):
		columns = 2
		rows = int(num_plots+1)/2
	else:
		rows = num_plots
		columns = 1
	xscale =  .4*dpi * 10*columns
	yscale = .4*dpi * 4 * rows
	#for i in xrange(0,int(f)+1,del_t*ndt):
	#	maxi(i,file_out,titles,nt,xmin,xmax,thresh_value)
	print file_out
	for j in xrange(num_plots):
		for j2 in xrange(len(file_out[j])):
			ind_fol = folder[j][j2]
			#max_min_thresh_data = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(maxi2)(i,[(file_out[j])[j2]],lineout[j]) for i in xrange(int((first_ind[j])[j2]),int((last_ind[j])[j2])+1,(frac[j])[j2]))
			max_min_thresh_data = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(maxi2)(i,[(file_out[j])[j2]],lineout[j]) for i in xrange(int(startn[ind_fol]),int(lastn[ind_fol])+1,int(ndump[ind_fol])))
			for i in xrange(0,len(max_min_thresh_data)):
				data_info = max_min_thresh_data[i]
				max_info = data_info[0]
				min_info = data_info[1]
				thresh_info = data_info[2]
				for jj2 in xrange(len(max_info)):
					if(max_info[jj2] > xmax[j]):
						xmax[j] = max_info[jj2]
					if(min_info[jj2] < xmin[j]):
						xmin[j] = min_info[jj2]
					if(thresh_info[jj2] > thresh_value[j]):
						thresh_value[j] = thresh_info[jj2]
	print xmax,xmin,thresh_value
	for jj in xrange(num_plots):
		if(max_min_overwrite[jj] == 'y'):
			xmax[jj] = max_over[jj]
			xmin[jj] = min_over[jj]
			if(log[jj] == 'y'):
				thresh_value[jj] = threshold_over[jj]
	Parallel(n_jobs=multiprocessing.cpu_count())(delayed(cock)(i,file_out,titles,1,xmin,xmax,thresh_value,log,colormap,lineout,index_info) for i in xrange(0,numx+1))
	#subprocess.call(['ffmpeg -framerate 10 -pattern_type glob -i plots/*.png -c:v libx264 -vf scale=-1:1080,format=yuv420p out_out2.mp4']) 
	
	plots = raw_input('Create movie? (y/n) : ')
	if(plots == 'y'):
		subprocess.call(["ffmpeg", "-framerate","10","-pattern_type","glob","-i",'plots/*.png','-c:v','libx264','-vf','scale=' + str(xscale) +':' + str(yscale)+' ,format=yuv420p','movie.mp4'])


def print_dirs(dirs):
	row_format ="{:>60}" +"{:>10}"
	lis =['Data']
	lis2 = ['Number','Data','Number']
	row_format2 = "{:>60}" +"{:>10}"+"{:>60}" +"{:>10}"
	print row_format2.format('Data',*lis2)
	print '------------------------------------------------------------------------------------------------------------------------------------------------'
	half = int((len(dirs)+1)/2)
	print len(dirs)
	pri = [(dirs)[:half],(dirs)[half:len(dirs)]]
	for j in xrange(half):
		if(j < half-1):
			print row_format2.format((pri[0])[j],str(j),(pri[1])[j],str(j+half))
		else:
			if(half * 2 == len(dirs)):
				print row_format2.format((pri[0])[j],str(j),(pri[1])[j],str(j+half))
			else:
				print row_format.format((pri[0])[j],str(j))

		print '------------------------------------------------------------------------------------------------------------------------------------------------'
	print ''

def cock(j,file_out,titles,nt,xmin,xmax,thresh_value,log,colormap,lineout,index_info):
	startn = index_info[0]
	ndump = index_info[1]
	lastn = index_info[2]
	folder = index_info[3]
	if(len(titles) > 3):
		columns = 2
		rows = int(len(titles)+1)/2
	else:
		rows = len(titles)
		columns = 1
	#print rows
	#print titles,file_out
	#print j 
	font_size = 13
	title_size = 20
	title_font_size = 14
	fig = plt.figure(1,figsize=(10*columns,4*rows))
	ti = nt*j
	print file_out
	ti = '%.2f' % ti
	plt.suptitle(r'$\/t\/ =\/' +ti+'\/[1/w_0]$',fontsize = title_size )
	for j2 in xrange(len(file_out)):
		ax = plt.subplot(rows,columns,j2+1)
		plt.title(titles[j2],fontsize = title_font_size)
		for jj2 in xrange(len(file_out[j2])):
			fol_num = folder[j2][jj2]
			filename = h5py.File(file_out[j2][jj2]+ str(str(int(1000000+j*ndump[fol_num]+startn[fol_num]))[1:])+'.h5', 'r')
			data=(filename[filename.keys()[1]])[:]
			#data2 =np.copy(data)
			#data2[data2 < 1e-3] = np.nan
			val_thres = thresh_value[j2]
			maxx = max(np.abs(xmin[j2]),xmax[j2])
			minn1 = min(np.abs(xmin[j2]),xmax[j2])
			if(log[j2] == 'y'):
				thresh = int(np.log10(val_thres))
				base = 10
				maxlog = int(np.log10(maxx))
				diff = maxlog -thresh
				if (diff > 4):
					del_base = 2
					diff = int(diff/2)*2
					ad_del = 1
					adjust_ind = 0
				else:
					del_base = 1
					ad_del = 0
					thresh = thresh +ad_del
					adjust_ind = 1
			#	if(xmax[j2] > 0):
			#		tick_locations=( [0.0] +[(maxx*base**x) for x in xrange(thresh,1,1)] )
			#		colmap = 'Reds'
			#		imAx = plt.imshow(data[ind:,:],aspect = 'auto',cmap = colmap,origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx, norm=matplotlib.colors.SymLogNorm( maxx*0.5))
			#	else:
			#		tick_locations=([-(maxx*base**x) for x in xrange(0,thresh-1,-1)] +[0.0] )
				#	colmap = 'Blues_r'
				ticklist = [-(maxx*base**x) for x in xrange(0,-diff+ad_del+adjust_ind,-del_base)]
				ticklist_reverse = ticklist[::-1]
				ticklist_reverse = [-1 * i for i in ticklist_reverse]	
				tick_locations=(ticklist +[0.0]+ticklist_reverse  )
				if(colormap[j2] == 'def'):
					imAx = plt.imshow(data[ind:,:],aspect = 'auto',origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx, norm=matplotlib.colors.SymLogNorm( 10**(thresh+1)))
				else:
					# imAx = plt.imshow(data[ind:,:],aspect = 'auto',origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx, norm=matplotlib.colors.SymLogNorm( 10**(thresh+1)),cmap=colormap[j2])
					imAx = plt.imshow(data[ind:,:],aspect = 'auto',origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx, norm=matplotlib.colors.SymLogNorm( 10**(thresh+1)),cmap=colormap[j2])

				divider = make_axes_locatable(ax)
				cax = divider.append_axes("right", size="2%", pad=0.05)
				plt.colorbar(imAx,cax=cax,ticks=tick_locations,format='%.1e')
			else:
				if(lineout[j2] == 'y'):
					plt.plot(data[300,:])
					axes = plt.gca()
					axes.set_ylim([xmin[j2],xmax[j2]])
				else:
					if(colormap[j2] == 'def'):
						imAx = plt.imshow(data[ind:,:],aspect = 'auto',origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx)
					else:
						imAx = plt.imshow(data[ind:,:],aspect = 'auto',origin='lower', interpolation='bilinear',vmin = -maxx,vmax = maxx,cmap=colormap[j2])
					divider = make_axes_locatable(ax)
					cax = divider.append_axes("right",size="2%",pad= 0.05)
					plt.colorbar(imAx,cax=cax,format='%.1e')

			filename.close()
	
	

	#plt.savefig('plots/charge'+str(j+1000000)[1:],dpi=400,bbox_inches = 'tight')
	plt.savefig('plots/charge'+str(j+1000000)[1:],dpi=dpi,bbox_inches = 'tight')
	plt.close()

def maxi(j,file_out,xmin,xmax,thresh_value):
	for j2 in xrange(len(file_out)):
		filename = h5py.File(file_out[j2]+ str(str(1000000+j)[1:])+'.h5', 'r')
		data=(filename[filename.keys[1]])[:]
		mx = np.max(data[ind:,:])
		mn = np.min(data[ind:,:])
		f = data[data != 0]
		if(f.size != 0):
			thresh = np.min(np.abs(f))
			if(thresh_value[j2] == 0):
				thresh_value[j2] = thresh
			else:
				thresh_value[j2] = np.min(thresh_value[j2],thresh)
		xmax[j2] = max(xmax[j2],mx)
		xmin[j2] = min(xmin[j2],mn)
		filename.close()

def maxi2(j,file_out,lineout):
	print file_out
	print j
	xmax_out = np.zeros(len(file_out))
	xmin_out = np.zeros(len(file_out))
	thresh_out = np.zeros(len(file_out))
	for j2 in xrange(len(file_out)):
		filename = h5py.File(file_out[j2] + str(str(1000000+j)[1:])+'.h5','r')
		data= (filename[filename.keys()[1]])[:]
		ln_out = lineout[j2]
		if(ln_out == 'y'):
			xmax_out[j2] = np.max(data[ind,:])
			xmin_out[j2] = np.min(data[ind,:])
			f = data[data != 0]
			if(f.size == 0):
					thresh_out[j2] =- np.inf
			else:
					thresh_out[j2] = np.min(np.abs(f))
		else:
			xmax_out[j2] = np.max(data[ind:,:])
			xmin_out[j2] = np.min(data[ind:,:])
			f = data[data != 0]
			if(f.size == 0):
				thresh_out[j2] =- np.inf
			else:
				thresh_out[j2] = np.min(np.abs(f))
	return np.array([xmax_out,xmin_out,thresh_out])
	

def fmt(x, pos):
	a, b = '{:.2e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)

def read_input_deck(file):
	file = open(file,'r')
	input_file = file.read()
	file.close()
	sim_params = read_simulation_params(input_file)
	num_plots = sim_params[0]
	data = []
	for j in xrange(1,num_plots+1):
		data_temp = read_subplot(input_file,j)
		data.append(data_temp)
	return [sim_params,data]
	
def read_simulation_params(string):
	ind = string.find('simulation')
	ind_start = string.find('{',ind)
	subplots_ind = string.find('subplots',ind_start)
	numsubplots = int(string[(string.find('=',subplots_ind)+1):string.find(',',subplots_ind)])

	nstart_ind = string.find('nstart',ind_start)
	nstart = int(string[(string.find('=',nstart_ind)+1):string.find(',',nstart_ind)])

	ndump_ind = string.find('ndump',ind_start)
	ndump = int(string[(string.find('=',ndump_ind)+1):string.find(',',ndump_ind)])

	nend_ind = string.find('nend',ind_start)
	nend = int(string[(string.find('=',nend_ind)+1):string.find(',',nend_ind)])

	return [numsubplots,nstart,ndump,nend]

def read_subplot(string,num):
	first_ind = 0
	for j in xrange(num):
		first_ind = string.find('data',first_ind)+4
	ind_start = string.find('{',first_ind)


	folders_ind = string.find('folders',ind_start)
	end_ind = string.find('\n',folders_ind)
	s2 = string[folders_ind:end_ind]
	quote_ind = s2.find('"')
	folders = []
	while(quote_ind != -1):
		quote_end = s2.find('"',quote_ind +1)
		folder_name = s2[(quote_ind+1):quote_end]
		folders.append(folder_name)
		quote_ind = s2.find('"',quote_end+1)

	title_ind = string.find('title',ind_start)
	quote_ind = string.find('"',title_ind)
	quote_end = string.find('"',quote_ind+1)
	title = string[quote_ind+1:quote_end]

	title_ind = string.find('logscale',ind_start)
	quote_ind = string.find('=',title_ind)
	quote_end = string.find(',',quote_ind+1)
	logscale = string[quote_ind+1:quote_end]
	logscale = logscale.strip() == 'True'


	title_ind = string.find('lineout',ind_start)
	quote_ind = string.find('=',title_ind)
	quote_end = string.find(',',quote_ind+1)
	lineout = string[quote_ind+1:quote_end]
	lineout = lineout.strip() == 'True'

	title_ind = string.find('min_max_overwrite',ind_start)
	quote_ind = string.find('=',title_ind)
	quote_end = string.find(',',quote_ind+1)
	overwrite = string[quote_ind+1:quote_end]
	overwrite = overwrite.strip() == 'True'

	if(not lineout):
		colormap_ind = string.find('colormap',ind_start)
		quote_ind = string.find('"',colormap_ind)
		quote_end = string.find('"',quote_ind+1)
		colormap = string[quote_ind+1:quote_end]
	else:
		colormap = 'def'
		


	if(overwrite):
		if(logscale):
			min_ind = string.find('minimum',ind_start)
			min_val = float(string[(string.find('=',min_ind)+1):string.find(',',min_ind)])
			max_ind = string.find('maximum',ind_start)
			max_val = float(string[(string.find('=',max_ind)+1):string.find(',',max_ind)])
			thresh_ind = string.find('threshold',ind_start)
			thresh_val = float(string[(string.find('=',thresh_ind)+1):string.find(',',thresh_ind)])
			axis_info = [min_val,max_val,thresh_val]
		else:
			min_ind = string.find('minimum',ind_start)
			min_val = float(string[(string.find('=',min_ind)+1):string.find(',',min_ind)])
			max_ind = string.find('maximum',ind_start)
			max_val = float(string[(string.find('=',max_ind)+1):string.find(',',max_ind)])
			axis_info = [min_val,max_val,0]
		return [folders,title,logscale,lineout,overwrite,axis_info,colormap]

	
		
			

	return [folders,title,logscale,lineout,overwrite,[-np.inf,np.inf,0],colormap]

	



if __name__ == "__main__":
	main()
