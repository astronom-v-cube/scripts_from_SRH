# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import scipy as SP;
import pylab as LAB;

def hhmm_format(t, pos):
  hh = (int)(t / 3600.);
  t -= hh*3600.;
  mm = (int)(t / 60.);
  return '%02d:%02d' % (hh,mm);

data = SP.io.readsav('C:/Users/Sergey/IDLWorkspace71/Default/flare20160306.sav');

#fig = LAB.figure();
#sp = fig.add_subplot(2,1,1);
#sp.xaxis.set_major_formatter(LAB.FuncFormatter(hhmm_format));

LAB.plot(data['map_time'],data['sunflux_h_lcp']);
LAB.plot(data['map_time'],data['sunflux_h_rcp']);

