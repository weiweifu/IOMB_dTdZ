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
    dTdZ = (v.data[i0,...]-v.data[ic,...])/vdiff*1e3
    dTdZ.shape = (1,) + dTdZ.shape
    
    slope = Variable(
                     time = np.asarray([0.]),
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
        
