from ILAMB.Confrontation import Confrontation
from ILAMB.Variable import Variable
import numpy as np

"""Notes

1) export PYTHONPATH=$PATH
2) add a time dimension to the obs, of size 1 where the bounds span
the entire time period of the GLODAP data.

     time = [(2000-1850)*365)]
     time_bounds = [[(1994-1850)*365),(2007-1850)*365)]]
     "days since 1850-01-01"

3) or use the keywords as below
4) to work with current ILAMB analysis, the returned 'slope' might
need to be temporal with a single time slice.

"""
def sc_fit(temp=None,modn='CESM2',X=None,Y=None,shift_xy=None,latbnd=slice(0,50),tsp='default',coef2d=None,y1=1.8,y2=4.5):
  from sklearn.linear_model import  LinearRegression

  tL={'default':slice(0,None),'1994s':slice(14,27),'1980s':slice(0,10),'1990s':slice(10,20),'2000s':slice(20,30)}
  def _dindx(v,S0=200,Sc=1000):
    # adep = v['lev_bnds'][:,1] - v['lev_bnds'][:,0]
    # find the index for 0 and 3000 m
    i0 = np.argmin(abs(v['lev']- S0))
    ic = np.argmin(abs(v['lev']- Sc))
    print('level ',i0,'-',ic,' in ',v['lev'][i0],'-',v['lev'][ic],'m')
    d1 = slice(i0,ic) # CMIP 0-3000m
    # return adep, d1
    return i0,ic
  def _getslope(X,Y3D):
    sh = Y3D.shape
    coef = copy.deepcopy(Y3D[0,:,:])
    print('the shape of Y3D is',sh)
    for j in range(0,sh[-2]):
      for i in range(0,sh[-1]):
        model = LinearRegression()
        Y = Y3D[:,j,i]
        output = model.fit(X[:,None],Y)
        coef[j,i] = output.coef_
    return coef

  if temp is not None:
    i0,ic = _dindx(temp)
    lev = temp['lev']
    X = lev[i0:ic+1]
    if temp['thetao'].ndim == 4:
      Y = temp['thetao'][tL[tsp],:,latbnd,:].mean(axis=(0,2,3))[i0:ic+1]
      Y3D = temp['thetao'][tL[tsp],i0:ic+1,:,:].mean(axis=(0))
    else:
      Y = temp['thetao'][:,latbnd,:].mean(axis=(1,2))[i0:ic+1]
  else:
    print('X and Y must be specified')
  print(X.shape,Y.shape)

  print('try to get coefs')
  if coef2d is not None:
    ccc = _getslope(X,Y3D)
  else:
    ccc = None

  model = LinearRegression()

  plt.close('all')
  fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(5,3))
  fig.tight_layout()
  if shift_xy is None:
    ax.plot(X,Y,'.')
    output = model.fit(X[:,None],Y)
  else:
    ax.plot(Y,X,'.')
    output = model.fit(Y[:,None],X)
  # print model
  print('Model intercept:',output.intercept_)
  print('Model slope:',output.coef_)
  if shift_xy is None:
    R2 = output.score(X[:,None],Y)
    print('R2:',R2)
    xfit = np.linspace(X.min(),X.max(),1000)
    yfit = output.predict(xfit[:,None])
    ax.plot(xfit,yfit,'r')
  else:
    R2 = output.score(X[:,None],Y)
    print('R2:',R2)
    xfit = np.linspace(Y.min(),Y.max(),1000)
    yfit = output.predict(xfit[:,None])
    ax.plot(xfit,yfit,'r')
  ax.xaxis.set_tick_params(labelsize=14)
  ax.yaxis.set_tick_params(labelsize=14)
  # ax.set_title(modn,fontsize=14)
  ax.text(np.median(xfit),np.median(yfit),'dT/dz='+'{:.4f}'.format(model.coef_[0])+'\nR2='+'{:.2f}'.format(R2),fontsize=14)
  ax.set_ylim([y1, y2])
  fig.suptitle(modn,fontsize=14)
  return fig,ccc

def GetSlope(v):
    # loop over unmasked cells and compute slope wrt depth
    def _dindx(v,S0=200,Sc=1000):
      # find the index for 200 and 1000 m
      dZ = np.asarray(getattr(v,'depth'))
      i0 = np.argmin(abs(dZ - S0))
      ic = np.argmin(abs(dZ - Sc))
      vdiff = dZ[ic] - dZ[i0]
      print('level {:2d}-{:2d} in {:6.2f}-{:6.2f}m time is {}'. 
                 format(i0,ic,dZ[i0],dZ[ic],v.time))
      return i0,ic,vdiff

    i0,ic,vdiff = _dindx(v)
    with np.errstate(under='ignore'):
        # dTdZ = (v.data[i0,...]-v.data[ic,...])/vdiff*1e3
        dTdZ = sc_fit(v)
    dTdZ.shape = (1,) + dTdZ.shape
    
    slope = Variable(
                     time = np.asarray([0.5]),
                     time_bnds = np.asarray([[0.,1.]]),
                     data = dTdZ,
                     lat       = v.lat,
                     lat_bnds  = v.lat_bnds,
                     lon       = v.lon,
                     lon_bnds  = v.lon_bnds,
                     unit = '1e-3 degC/m'
                     )
    return slope 

class ConfdTdZ(Confrontation):

    def stageData(self,m):

        y0 = float(self.keywords.get("y0",2000.)) # [yr] beginning year to include in analysis
        yf = float(self.keywords.get("yf",2001.)) # [yr] end year to include in analysis
        
        obs = Variable(filename       = self.source,
                       variable_name  = self.variable,
                       alternate_vars = self.alternate_vars,
                       t0 = None if len(self.study_limits) != 2 else self.study_limits[0],
                       tf = None if len(self.study_limits) != 2 else self.study_limits[1])
        #if obs.time is None: raise il.NotTemporalVariable()
        self.pruneRegions(obs)

        mod = m.extractTimeSeries(self.variable,
                                  alt_vars     = self.alternate_vars,
                                  expression   = self.derived,
                                  initial_time = (y0-1850)*365,
                                  final_time   = (yf-1850)*365,
                                  lats         = None if obs.spatial else obs.lat,
                                  lons         = None if obs.spatial else obs.lon)
        
        # remove the time dimension
        obs = obs.integrateInTime(mean=True)
        mod = mod.integrateInTime(mean=True)

        # get slope
        obs_slope = GetSlope(obs)
        mod_slope = GetSlope(mod)

        return obs_slope,mod_slope
        
