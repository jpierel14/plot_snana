from __future__ import print_function
import numpy as np
import matplotlib as mpl
mpl.use('tkAgg')
import matplotlib.pyplot as plt
import os,glob,math,sys
from optparse import OptionParser
from scipy.interpolate import interp1d
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
from copy import copy

__band_order__=np.append(['u','b','g','r','i','z','y','j','h','k'],
	 [x.upper() for x in ['u','b','g','r','i','z','y','j','h','k']])

class snana_dat_plot(QtGui.QWidget):
	def __init__(self,cid_list,all_light_curves=None,all_spec_data=None,all_fits=None,all_peaks=None,bin_size=0.0):
		super(snana_dat_plot, self).__init__()
		
		
		self.bin_size=bin_size
		self.all_fits=all_fits
		self.all_peaks=all_peaks
		self.CIDs=cid_list
		self.CID=cid_list[0]
		if all_light_curves is not None:
			self.all_light_curves=all_light_curves
			self.light_curve_data=self.all_light_curves[self.CID]
		else:
			self.all_light_curves=None
			self.light_curve_data=None
		if all_spec_data is not None:
			self.all_spec_data=all_spec_data
			self.spec_data=self.all_spec_data[self.CID]
		else:
			self.all_spec_data=all_spec_data
			self.spec_data=None
		if all_fits is not None:
			self.fit=all_fits[self.CID]
		else:
			self.fit=None
		if all_peaks is not None:
			self.peak=all_peaks[self.CID]
		else:
			self.peak=None
		self.widgets=[]
		self.resize(1800,1800)
		self.init_ui()
		self.qt_connections()
		self.CID=cid_list[0]
		self.updateplot()

	def init_ui(self):
		
		layout = QtGui.QGridLayout()
		self.setLayout(layout)
		if self.all_light_curves is not None:
			self.setWindowTitle('SNANA LC Plot')
			for i in range(len(np.unique(self.light_curve_data['filter']))):
				self.widgets.append([pg.PlotWidget(),pg.PlotWidget()])
				
				layout.addWidget(self.widgets[i][0],i,0)
				layout.addWidget(self.widgets[i][1],i,1)

		elif self.all_spec_data is not None:
			self.setWindowTitle('SNANA SPEC Plot')
			for i in range(len(np.unique(self.spec_data['tobs']))):
				self.widgets.append([pg.PlotWidget(),pg.PlotWidget()])
				
				layout.addWidget(self.widgets[i][0],i,0)
				layout.addWidget(self.widgets[i][1],i,1)
		else:
			raise RuntimeError("No data to plot.")
		
		
		self.buttons={}
		j=0
		for cid in self.CIDs:
			temp_button=QtGui.QPushButton(cid)
			layout.addWidget(temp_button,i+1,j)
			self.buttons[cid]=temp_button
			j+=1
		#reset=QtGui.QPushButton('Reset')
		#layout.addWidget(reset,i+1,j)
		#self.reset=reset
		#self.setGeometry(10, 10, 1000, 600)
		self.show()

	def qt_connections(self):
		self.button_func_dict={}
		for cid in self.CIDs:
			self.button_func_dict[cid]=self.on_cid_button_clicked_generator(cid)
		for cid in self.CIDs:
			self.buttons[cid].clicked.connect(self.button_func_dict[cid])

	#def reset_plot(self):
	#self.ViewBox.clear()
	#	self.updateplot()

	def updateplot(self):
		if self.all_light_curves is not None and self.all_spec_data is None:
			for i in range(len(np.unique(self.light_curve_data['filter']))):
				for j in range(2):
					self.widgets[i][j].clear()
			self.plot_lc()
		else:
			for i in range(len(np.unique(self.spec_data['tobs']))):
				for j in range(2):
					self.widgets[i][j].clear()
			self.plot_spec()


	def on_cid_button_clicked_generator(self,cid):
		def button_func():
			print ("Switched to SN%s"%cid)
			self.CID=cid
			self.light_curve_data = self.all_light_curves[self.CID] if self.all_light_curves is not None else None
			self.spec_data=self.all_spec_data[self.CID] if self.all_spec_data is not None else None
			self.updateplot()

			
		return(button_func)
		

	def plot_lc(self):
		
		for i in range(len(np.unique(self.light_curve_data['filter']))):

			band=np.append([x for x in __band_order__ if x in np.unique(self.light_curve_data['filter'])],
							[x for x in np.unique(self.light_curve_data['filter']) if x not in __band_order__])[i]
			temp_sn={k:self.light_curve_data[k][np.where(self.light_curve_data['filter']==band)[0]] for k in self.light_curve_data.keys()}
			chi2=np.mean(temp_sn['chi2'])
			if chi2>0:
				lab=r'%s: $\chi^2$=%.1f'%(band,np.mean(temp_sn['chi2']))
				leg_size=12
			else:
				lab=band
				leg_size=16

			#p = self.plotcurve.addPlot(title='SN%s'%self.CID)
		
			scat=pg.ScatterPlotItem(temp_sn['time'],temp_sn['flux'],symbol='o', 
				pen=None,symbolPen=None, symbolSize=20, symbolBrush='w',title=band+' Band')
			#ax=pg.LabelItem(lab['bottom'])
			#ax.setParentItem(scat)
			
			err = pg.ErrorBarItem(x=temp_sn['time'], y=temp_sn['flux'], height=temp_sn['fluxerr'], beam=0.25, pen={'color':'w', 'width':.4})
			#leg = pg.LegendItem(size=(10,10),offset=(10,10))
			#leg.setParentItem(scat)
			#leg.addItem(scat,band)
			leg=pg.TextItem(band+' Band','w')
			self.widgets[i][0].addItem(err)
			self.widgets[i][0].addItem(scat)
			self.widgets[i][0].addItem(leg)
			ax=self.widgets[i][0].getAxis('left')
			ax.setLabel('Flux')

			#self.widgets[i][0].addItem(ax)
			if i ==0:
				lr = pg.LinearRegionItem([temp_sn['time'].min(),temp_sn['time'].max()])
				lr.setZValue(-10)
				self.widgets[i][0].addItem(lr)

			
			#self.plotcurve2.setData(self.light_curve_data['time'],self.light_curve_data['flux'],symbol='o',pen=None, symbolPen=None, symbolSize=20, symbolBrush='w')
			err2 = pg.ErrorBarItem(x=temp_sn['time'], y=temp_sn['flux'], height=temp_sn['fluxerr'], beam=0.25, pen={'color':'w', 'width':.4})
			scat2=pg.ScatterPlotItem(temp_sn['time'],temp_sn['flux'],symbol='o', pen=None,symbolPen=None, symbolSize=20, symbolBrush='w',name=band)
			leg2=pg.TextItem(band+' Band','w')
			self.widgets[i][1].addItem(err2)
			self.widgets[i][1].addItem(scat2)
			self.widgets[i][1].addItem(leg2)
			ax=self.widgets[i][1].getAxis('left')
			ax.setLabel('Flux')
			def updatePlot():
				for j in range(len(np.unique(self.light_curve_data['filter']))):
					self.widgets[j][1].setXRange(*lr.getRegion(), padding=0)
			def updateRegion():
				lr.setRegion(self.widgets[0][1].getViewBox().viewRange()[0])
			lr.sigRegionChanged.connect(updatePlot)
			self.widgets[i][1].sigXRangeChanged.connect(updateRegion)
			if len(self.fit)>0:
				fit_time=np.arange(temp_sn['time'][0],temp_sn['time'][-1],1)
				self.plotcurve.setData(fit_time,self.fit[band](fit_time),pen='r',width=3)
			
		ax=self.widgets[i][0].getAxis('bottom')
		ax.setLabel('Time-%.2f (Rest Frame Days)'%self.peak)
		ax=self.widgets[i][1].getAxis('bottom')
		ax.setLabel('Time-%.2f (Rest Frame Days)'%self.peak)
		
	def plot_spec(self):
		for j in range(len(np.unique(self.spec_data['tobs']))):

			temp_sn=np.where(self.spec_data['tobs']==np.unique(self.spec_data['tobs'])[j])[0]
			if self.bin_size!=0:
				
				binned_wave=[]
				binned_flux=[]
				binned_fluxerr=[]
				bins=np.trunc(self.spec_data['wave'][temp_sn]/bin_size)
				for i in np.unique(bins):
					binned_wave=np.append(binned_wave,np.mean(self.spec_data['wave'][temp_sn][bins==i]))
					binned_flux=np.append(binned_flux,np.mean(self.spec_data['flux'][temp_sn][bins==i]))
					binned_fluxerr=np.append(binned_fluxerr,np.mean(self.spec_data['fluxerr'][temp_sn][bins==i]))
			else:
				binned_wave=self.spec_data['wave'][temp_sn]
				binned_flux=self.spec_data['flux'][temp_sn]
				binned_fluxerr=self.spec_data['fluxerr'][temp_sn]
			binned_fluxerr/=np.max(binned_flux)
			binned_flux/=np.max(binned_flux)
			if j ==0:
				lr = pg.LinearRegionItem([binned_wave.min(),binned_wave.max()])
				lr.setZValue(-10)
				self.widgets[j][0].addItem(lr)
			line=pg.PlotCurveItem(binned_wave,binned_flux,pen='w')
			
			upper=pg.PlotCurveItem(binned_wave,binned_flux+binned_fluxerr,pen='r')
			lower=pg.PlotCurveItem(binned_wave,binned_flux-binned_fluxerr,pen='r')
			err=pg.FillBetweenItem(lower,upper,brush='r')#(255,0,0,.5))

			leg=pg.TextItem('Tobs: %.2f'%np.unique(self.spec_data['tobs'])[j],'w')
			leg2=pg.TextItem('Tobs: %.2f'%np.unique(self.spec_data['tobs'])[j],'w')
			self.widgets[j][0].addItem(err)
			self.widgets[j][0].addItem(line)
			
			self.widgets[j][0].addItem(leg)
			
			ax=self.widgets[j][0].getAxis('left')
			ax.setLabel('Normalized Flux')

			line2=pg.PlotCurveItem(binned_wave,binned_flux,pen='w')
			upper2=pg.PlotCurveItem(binned_wave,binned_flux+binned_fluxerr,pen='r')
			lower2=pg.PlotCurveItem(binned_wave,binned_flux-binned_fluxerr,pen='r')
			err2=pg.FillBetweenItem(lower,upper,brush='r')#(255,0,0,.5))
			self.widgets[j][1].addItem(err2)
			self.widgets[j][1].addItem(line2)
			self.widgets[j][1].addItem(leg2)
			
			ax=self.widgets[j][1].getAxis('left')
			ax.setLabel('Normalized Flux')
			def updatePlot():
				for k in range(len(np.unique(self.spec_data['tobs']))):
					self.widgets[k][1].setXRange(*lr.getRegion(), padding=0)
			def updateRegion():
				lr.setRegion(self.widgets[0][1].getViewBox().viewRange()[0])
			lr.sigRegionChanged.connect(updatePlot)
			self.widgets[j][1].sigXRangeChanged.connect(updateRegion)
			
			
		ax=self.widgets[j][0].getAxis('bottom')
		ax.setLabel('Wavelength (Angstrom)')
		ax=self.widgets[j][1].getAxis('bottom')
		ax.setLabel('Wavelength (Angstrom)')
	
	



def read_spec(cid,base_name):
	names=['wave','flux','fluxerr','tobs']
	id_to_obs=dict([])
	with open(base_name+".SPECLIST.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		else:
			id_to_obs[int(temp[varnames.index('ID')])]=float(temp[varnames.index('TOBS')])
	sn={k:[] for k in names}

	with open(base_name+".SPECPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		
		
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		elif len(temp)>0 and b'OBS:' in temp and\
			 str(temp[varnames.index('CID')].decode('utf-8'))in cid:
			sn['wave'].append((float(temp[varnames.index('LAMMAX')])+float(temp[varnames.index('LAMMIN')]))/2.)
			sn['flux'].append(float(temp[varnames.index('FLAM')]))
			sn['fluxerr'].append(float(temp[varnames.index('FLAMERR')]))
			sn['tobs'].append(id_to_obs[int(temp[varnames.index('ID')])])
	sn={k:np.array(sn[k]) for k in sn.keys()}
	return(sn)
def read_lc(cid,base_name):
	names=['time','flux','fluxerr','filter','chi2']
	peak=None
	sn={k:[] for k in names} 
	fit={k:[] for k in ['time','flux','filter']}
	with open(base_name+".LCPLOT.TEXT",'rb') as f:
		dat=f.readlines()
	for line in dat:
		temp=line.split()
		if len(temp)>0 and b'VARNAMES:' in temp:
			varnames=[str(x.decode('utf-8')) for x in temp]
		elif len(temp)>0 and b'OBS:' in temp and str(temp[varnames.index('CID')].decode('utf-8')) in cid:
			if int(temp[varnames.index('DATAFLAG')])==1:
				if peak is None:
					peak=float(temp[varnames.index('MJD')])-float(temp[varnames.index('Tobs')])
				sn['time'].append(float(temp[varnames.index('Tobs')]))
				sn['flux'].append(float(temp[varnames.index('FLUXCAL')]))
				sn['fluxerr'].append(float(temp[varnames.index('FLUXCAL_ERR')]))
				sn['filter'].append(str(temp[varnames.index('BAND')].decode('utf-8')))
				sn['chi2'].append(float(temp[varnames.index('CHI2')]))
			elif int(temp[varnames.index('DATAFLAG')])==0:
				fit['time'].append(float(temp[varnames.index('Tobs')]))
				fit['flux'].append(float(temp[varnames.index('FLUXCAL')]))
				fit['filter'].append(str(temp[varnames.index('BAND')].decode('utf-8')))
	
	sn={k:np.array(sn[k]) for k in sn.keys()}
	fit={k:np.array(fit[k]) for k in fit.keys()}
	if len(fit['filter'])>0:
		fits={k:interp1d(fit['time'][fit['filter']==k],
					 fit['flux'][fit['filter']==k]) for k in np.unique(fit['filter'])}
	else:
		fits=[]
	return(sn,fits,peak)

def get_spec_data(cid,base_name):
	sn=read_spec(cid,base_name)
	return(sn)
	



def get_lc_data(cid,base_name):
	sn,fits,peak=read_lc(cid,base_name)
	return(sn,fits,peak)
	

def plot_cmd(genversion,cid_list):
	if os.path.splitext(genversion)[1]=='.NML':
		plotter='salt2'
	else:
		plotter='normal'
	rand=str(np.random.randint(10000,100000))
	cmd="snana.exe NOFILE VERSION_PHOTOMETRY "+genversion+\
		" SNCID_LIST "+cid_list+\
		" CUTWIN_CID 0 0 SNTABLE_LIST 'SNANA(text:key) LCPLOT(text:key) SPECPLOT(text:key)' TEXTFILE_PREFIX 'OUT_TEMP_"+rand+\
		"' > OUT_TEMP_"+rand+".LOG"
	os.system(cmd)
	return(plotter,'OUT_TEMP_'+rand)

def main():
	parser = OptionParser()
	parser.add_option("--spec",action="store_true",dest="spec",default=False)
	parser.add_option("--lc",action="store_true",dest="lc",default=False)
	parser.add_option("-i",action="store",type="string",dest="CID",default="None")
	parser.add_option("-b",action="store",type="float",dest='bin_size',default=0)
	parser.add_option("-g",action="store",type='string',dest='genversion',default=None)
	(options,args)=parser.parse_args()
	if options.CID=="None":
		raise RuntimeError("Need to define CID")
	if options.genversion is None:
		raise RuntimeError("Need to define genversion")
	

	plotter_choice,options.base_name=plot_cmd(options.genversion,options.CID)
	options.CID=options.CID.split(',')
	#QtGui.QApplication.setGraphicsSystem('raster')
	app = QtGui.QApplication([])
	#mw = QtGui.QMainWindow()
	#mw.resize(800,800)

	#win = pg.GraphicsWindow(title="Basic plotting examples")
	#win.resize(1000,600)
	#win.setWindowTitle('pyqtgraph example: Plotting')

	# Enable antialiasing for prettier plots
	pg.setConfigOptions(antialias=True)


	
	
	if options.spec:
		all_spec={}
		for cid in options.CID:
			all_spec[cid]=get_spec_data(cid,options.base_name)
		temp=snana_dat_plot(options.CID,all_spec_data=all_spec,bin_size=options.bin_size)
		
	elif options.lc:
		all_lc={}		
		all_fits={}
		all_peak={}
		for cid in options.CID:
			all_lc[cid],all_fits[cid],all_peak[cid]=get_lc_data(cid,options.base_name)
		temp=snana_dat_plot(options.CID,all_light_curves=all_lc,all_fits=all_fits,all_peaks=all_peak)
		
	else:
		all_lc={}
		all_spec={}
		all_fits={}
		all_peak={}
		for cid in options.CID:
			all_lc[cid],all_fits[cid],all_peak[cid]=get_lc_data(cid,options.base_name)
			all_spec[cid]=get_spec_data(cid,options.base_name)
		temp=snana_dat_plot(options.CID,all_spec_data=all_spec,bin_size=options.bin_size)
		temp2=snana_dat_plot(options.CID,all_light_curves=all_lc,all_fits=all_fits,all_peaks=all_peak)
	try:
		QtGui.QApplication.instance().exec_()
		for x in glob.glob(options.base_name+'*'):
			os.remove(x)
	except RuntimeError:
		for x in glob.glob(options.base_name+'*'):
			os.remove(x)

if __name__=='__main__':
	main()
