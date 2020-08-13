from ILAMB.Confrontation import Confrontation
from ILAMB.Variable import Variable
import numpy as np
import ntpath

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
class ConfCant(Confrontation):

    def stageData(self,m):

        y0 = float(self.keywords.get("y0",1994.)) # [yr] beginning year to include in analysis
        yf = float(self.keywords.get("yf",2007.)) # [yr] end year to include in analysis
        
        obs = Variable(filename       = self.source,
                       variable_name  = self.variable,
                       alternate_vars = self.alternate_vars)
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

        return obs,mod
        
