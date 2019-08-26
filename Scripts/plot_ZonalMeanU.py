"""
Script plots the zonal mean zonal wind response across all latitudes.
Statistical test uses the FDR method with alpha_FDR=0.05

Notes
-----
    Author : Zachary Labe
    Date   : 26 August 2019
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
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
print('\n' '----Plotting Zonal Mean U - %s----' % titletime)

### Alott time series (100 ensemble members)
year1 = 1901
year2 = 2000
years = np.arange(year1,year2+1,1)

###############################################################################
###############################################################################
###############################################################################
### Call arguments
letters = list(string.ascii_lowercase)

### Define directories
directorydata = '/seley/zlabe/simu/'
directoryfigure = '/home/zlabe/Desktop/SSTVARI/ZonalMeanU/'

######################
def readDataPeriods(varnames,simulations,period):
    ### Call function for 4d variable data
    lat,lon,lev,varfuture = MO.readExperiAll(varnames,simulations[0],'zonmean')
    lat,lon,lev,varpast = MO.readExperiAll(varnames,simulations[1],'zonmean')
    
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
        print(varpast.shape)
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
    anom = varfuturem - varpastm
    
    ### Calculate ensemble mean
    anommean = np.nanmean(anom,axis=0)
    
    ### Calculate climatologies
    climo = np.nanmean(varpastm,axis=0)
    
    ### Calculate significance for each month (pick method)
    pruns = UT.calc_FDR_ttest(varfuturem[:,:,:],varpastm[:,:,:],0.05) #FDR
    
    return anommean,climo,pruns,lat,lon,lev

### Read in data
anomy,climoy,pvaly,lat,lon,lev = readDataPeriods('U',['SST_Fu','SST_Cu'],'Annual')
anomwi,climowi,pvalwi,lat,lon,lev = readDataPeriods('U',['SST_Fu','SST_Cu'],'JFM')
anomsu,climosu,pvalsu,lat,lon,lev = readDataPeriods('U',['SST_Fu','SST_Cu'],'JJA')

data = [anomy,anomwi,anomsu]
clim = [climoy,climowi,climosu]
pp = [pvaly,pvalwi,pvalsu]

###########################################################################
###########################################################################
###########################################################################
#### Plot U
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Set limits for contours and colorbars
limit = np.arange(-4,4.1,0.25)
barlim = np.arange(-4,5,1)
zscale = np.array([1000,700,500,300,200,
                    100,50,30,10])
latq,levq = np.meshgrid(lat,lev)

fig = plt.figure(figsize=(6,8))
for i in range(len(data)):
    anom = data[i]
    climo = clim[i]
    pval = pp[i]
    
    ax1 = plt.subplot(3,1,i+1)
    
    if i == 0:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelbottom=False)    
        ax1.yaxis.set_ticks_position('left')
    elif i == 1:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelbottom=False)    
    elif i == 2:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
    elif i == 3:
        ax1.spines['top'].set_color('dimgrey')
        ax1.spines['right'].set_color('dimgrey')
        ax1.spines['bottom'].set_color('dimgrey')
        ax1.spines['left'].set_color('dimgrey')
        ax1.spines['left'].set_linewidth(2)
        ax1.spines['bottom'].set_linewidth(2)
        ax1.spines['right'].set_linewidth(2)
        ax1.spines['top'].set_linewidth(2)
        ax1.tick_params(axis='y',direction='out',which='major',pad=3,
                        width=0,color='dimgrey',labelleft=False)
        ax1.tick_params(axis='x',direction='out',which='major',pad=3,
                        width=2,color='dimgrey')    
        ax1.xaxis.set_ticks_position('bottom')
    
    ### Set limits    
    cs = plt.contourf(lat,lev,anom,limit,extend='both')
    cs1 = plt.contour(lat,lev,climo,np.arange(0,81,5),
                      linewidths=1,colors='k') 
    cs2 = plt.contourf(lat,lev,pval,
                       colors='None',hatches=['////'])     
    
    cs.set_cmap(cmocean.cm.balance)
    
    plt.yscale('log',nonposy='clip')
    plt.xlim([-90,90])
    plt.ylim([1000,10])
    plt.xticks(np.arange(-90,96,15),map(str,np.arange(-90,91,15)),fontsize=8)
    plt.yticks(zscale,map(str,zscale),ha='right',fontsize=6)
    plt.minorticks_off()
    plt.gca().invert_yaxis()
    
    plt.ylabel(r'\textbf{Latitude [$\bf{^{\circ}}$N]}',fontsize=6,
                     color='dimgray')    
    ax1.annotate(r'\textbf{[%s]}' % letters[i],xy=(0,0),
            xytext=(0.01,0.84),xycoords='axes fraction',
            color='k',fontsize=11)

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=False)
cbar.set_label(r'\textbf{U [m/s]}',fontsize=11,color='dimgray')
cbar.set_ticks(barlim)
cbar.set_ticklabels(list(map(str,barlim))) 
cbar.ax.tick_params(axis='x', size=.01)
cbar.outline.set_edgecolor('dimgrey')
cbar.ax.tick_params(labelsize=6)

plt.tight_layout()
plt.subplots_adjust(bottom=0.2,hspace=0.06,wspace=0.08,right=0.95)

plt.savefig(directoryfigure + 'U_zonalmean.png',dpi=300)