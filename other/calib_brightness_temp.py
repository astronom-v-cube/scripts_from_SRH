from functools import lru_cache
import numpy as np

class mw_spectrum():
    def __init__(self, filename):
        self._filename = filename
    
    @lru_cache()
    def _read_file(self):
        print("Open file and read to numpy array")
        return np.genfromtxt(self._filename, delimiter=' ').T

    def get_btemp(self, frequency):
        data = self._read_file()
        idx = np.argmin(np.abs(data[0]- frequency))
        return data[1][idx]

ratan = mw_spectrum('grechnev_ratan_fit.dat')

# put frequency in GHz
print(ratan.get_btemp(2.802))
print(ratan.get_btemp(5.7))
print(ratan.get_btemp(17))
