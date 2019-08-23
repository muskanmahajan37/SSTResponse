"""
Plot map of PAMIP data for data comparing different simulations. 
Statistical test uses the FDR method with alpha_FDR=0.05

Notes
-----
    Author : Zachary Labe
    Date   : 23 August 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_MonthlyData as MO
import calc_Utilities as UT
import cmocean
import string

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print('\n' '----Plotting Monthly Maps - %s----' % titletime)

### Alott time series (100 ensemble members)
year1 = 1901
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
varnames = ['U10','Z50','U200','U700','Z500','SLP','P','T2M','RNET','SST']
letters = list(string.ascii_lowercase)
readallinfo = True
period = 'DJF'

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/SSTVARI/GlobalMaps/FuCu/%s/' % period

######################
def readDataPeriods(varnames,simulations,period):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,simulations[0],'surface')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,simulations[1],'surface')
    
    ### Create 2d array of latitude and longitude
    lon2,lat2 = np.meshgrid(lon,lat)
    
    ### Remove missing data [ensembles,months,lat,lon]
    varfuture[np.where(varfuture <= -1e10)] = np.nan
    varpast[np.where(varpast <= -1e10)] = np.nan
    
    ### Calculate over DJF
    if period == 'JJA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,5:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,5:8:,:,:],axis=1)
    elif period == 'JAS':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:9,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:9:,:,:],axis=1)
    elif period == 'July':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:7,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:7:,:,:],axis=1)
    elif period == 'August':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,7:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,7:8:,:,:],axis=1)
    elif period == 'JA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,6:8,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,6:8:,:,:],axis=1)
    elif period == 'S':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,8:9,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,8:9:,:,:],axis=1)
    elif period == 'AMJ':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,3:6,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,3:6:,:,:],axis=1)
    elif period == 'OND':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,-3:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,-3:,:,:],axis=1)
    elif period == 'DJF':
        print('Calculating over %s months!' % period)
        runs = [varfuture,varpast]
        var_mo = np.empty((2,varpast.shape[0]-1,varpast.shape[2],varpast.shape[3]))
        for i in range(len(runs)):
            var_mo[i,:,:,:] = UT.calcDecJanFeb(runs[i],runs[i],lat,lon,'surface',1) 
        varfuturem = var_mo[0]
        varpastm = var_mo[1]
    elif period == 'JFM':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:3,:,:],axis=1)
    elif period == 'JF':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:2,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:2,:,:],axis=1)
    elif period == 'FMA':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:4,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:4,:,:],axis=1)
    elif period == 'FM':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:3,:,:],axis=1)
    elif period == 'J':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,0:1,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,0:1,:,:],axis=1)
    elif period == 'F':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,1:2,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,1:2,:,:],axis=1)
    elif period == 'M':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,2:3,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,2:3,:,:],axis=1)
    elif period == 'Annual':
        print('Calculating over %s months!' % period)
        varfuturem = np.nanmean(varfuture[:,:,:,:],axis=1)
        varpastm = np.nanmean(varpast[:,:,:,:],axis=1)
    else:
        print(ValueError('Selected wrong month period!'))
    
    ### Calculate anomalies
    anompi = varfuturem - varpastm
    
    ### Calculate ensemble mean
    anompim = np.nanmean(anompi,axis=0)
    zdiffruns = anompim
    
    ### Calculate climatologies
    zclimo = np.nanmean(varpastm,axis=0)
    
    ### Calculate significance for each month (pick method)
    pruns = UT.calc_FDR_ttest(varfuturem[:,:,:],varpastm[:,:,:],0.05) #FDR
#    pruns = UT.calc_indttest(varfuturem[:,:,:],varpastm[:,:,:])[1] #ttest
    
    return zdiffruns,zclimo,pruns,lat,lon

###########################################################################
###########################################################################
###########################################################################
### Read in data
if readallinfo == True:
    diff = np.empty((len(varnames),96,144))
    climo = np.empty((diff.shape))
    pp = np.empty((diff.shape))
    lat = np.empty((len(varnames),96))
    lon = np.empty((len(varnames),144))
    for v in range(len(varnames)):
        diff[v],climo[v],pp[v],lat[v],lon[v] = readDataPeriods(
                                                        varnames[v],
                                                         ['SST_Fu','SST_Cu'],
                                                         period)
 
##########################################################################
##########################################################################
##########################################################################
### Plot settings
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

for i in range(len(varnames)):
    print('Completed: Plot for %s!' % (varnames[i]))
    variables = diff[i,:,:]
    climos = climo[i,:,:] 
    pvalues = pp[i,:,:]
    
    ### Plot settings
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
    
    ### Set limits for contours and colorbars
    if varnames[i] == 'T2M':
        limit = np.arange(-10,10.1,0.5)
        barlim = np.arange(-10,11,10)
    if varnames[i] == 'SST':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnames[i] == 'SLP':
        limit = np.arange(-4,4.1,0.25)
        barlim = np.arange(-4,5,4)
    elif varnames[i] == 'Z500':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50) 
    elif varnames[i] == 'Z50':
        limit = np.arange(-100,100.1,10)
        barlim = np.arange(-100,101,100) 
    elif varnames[i]=='U10':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnames[i]=='U200':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnames[i]=='U700':
        limit = np.arange(-3,3.1,0.25)
        barlim = np.arange(-3,4,3)
    elif varnames[i] == 'SWE':
        limit = np.arange(-25,25.1,1)
        barlim = np.arange(-25,26,25)
    elif varnames[i] == 'P':
        limit = np.arange(-2,2.1,0.05)
        barlim = np.arange(-2,3,2) 
    elif varnames[i] == 'THICK':
        limit = np.arange(-60,60.1,3)
        barlim = np.arange(-60,61,60)
    elif varnames[i] == 'EGR':
        limit = np.arange(-0.2,0.21,0.02)
        barlim = np.arange(-0.2,0.3,0.2)
    elif varnames[i] == 'RNET':
        limit = np.arange(-50,50.1,2)
        barlim = np.arange(-50,51,50)
        
    ### Meshgrid lat and lon    
    lonq,latq = np.meshgrid(lon,lat)
    
    fig = plt.figure()
    ax1 = plt.subplot(111)

    m = Basemap(projection='moll',lon_0=0,resolution='l') 
    
    varn, lons_cyclic = addcyclic(variables, lon[i])
    varn, lons_cyclic = shiftgrid(180., varn, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lat[i])
    x, y = m(lon2d, lat2d)
    
    pvarn,lons_cyclic = addcyclic(pvalues, lon[i])
    pvarn,lons_cyclic = shiftgrid(180.,pvarn,lons_cyclic,start=False)
    climoq,lons_cyclic = addcyclic(climos, lon[i])
    climoq,lons_cyclic = shiftgrid(180.,climoq,lons_cyclic,start=False)
              
    circle = m.drawmapboundary(fill_color='white',
                               color='dimgrey',linewidth=0.7)
    circle.set_clip_on(False)
    
    if varnames[i] == 'RNET':
        varn = varn * -1. # change sign for upward fluxes as positive
    
    ### Plot contours
    cs = m.contourf(x,y,varn,limit,extend='both')
    cs1 = m.contourf(x,y,pvarn,colors='None',hatches=['.....'])
    if varnames[v] == 'Z50': # the interval is 250 m 
        cs2 = m.contour(x,y,climoq,np.arange(21900,23500,250),
                        colors='k',linewidths=1.1,zorder=10)
        
    ### Set map geography
    if varnames[i] == 'RNET' or varnames[i] == 'SST':
        m.drawcoastlines(color='darkgrey',linewidth=0.15)
        m.fillcontinents(color='dimgrey')
    else:
        m.drawcoastlines(color='dimgrey',linewidth=0.5)
   
    ### Set colormap
    cs.set_cmap(cmocean.cm.balance)
    
    plt.title(r'\textbf{%s}' % varnames[i],color='dimgrey',fontsize=30)
    
    cbar_ax = fig.add_axes([0.312,0.07,0.4,0.022])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='max',extendfrac=0.07,drawedges=True)
        
    if varnames[i] == 'T2M' or varnames[i] == 'SST':
        cbar.set_label(r'\textbf{$^\circ$C}',fontsize=7,color='dimgray')  
    elif varnames[i] == 'Z500':
        cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray')  
    elif varnames[i] == 'Z50':
        cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray')  
    elif varnames[i] == 'SLP':
        cbar.set_label(r'\textbf{hPa}',fontsize=7,color='dimgray')  
    elif varnames[i] == 'U10' or varnames[i] == 'U200' or varnames[i] == 'U500' or varnames[v] == 'U700':
        cbar.set_label(r'\textbf{m/s}',fontsize=7,color='dimgray')  
    elif varnames[i] == 'SWE':
        cbar.set_label(r'\textbf{mm}',fontsize=7,color='dimgray')
    elif varnames[i] == 'P':
        cbar.set_label(r'\textbf{mm/day}',fontsize=7,color='dimgray') 
    elif varnames[i] == 'THICK':
        cbar.set_label(r'\textbf{m}',fontsize=7,color='dimgray') 
    elif varnames[i] == 'EGR':
        cbar.set_label(r'\textbf{1/day}',fontsize=7,color='dimgray')
    elif varnames[i] == 'RNET':
        cbar.set_label(r'\textbf{W/m$\bf{^{2}}$',fontsize=7,color='dimgray')
            
    cbar.set_ticks(barlim)
    cbar.set_ticklabels(list(map(str,barlim))) 
    cbar.ax.tick_params(axis='x', size=.01)
    cbar.outline.set_edgecolor('dimgrey')
    cbar.outline.set_linewidth(0.5)
    cbar.ax.tick_params(labelsize=5)

    plt.tight_layout()
    
    plt.savefig(directoryfigure + '%s_%s_FDR.png' % (varnames[i],period),
                dpi=300)
    print('Completed: Script done!')


        