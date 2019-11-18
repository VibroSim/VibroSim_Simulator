# enter_frequency.py: Convenience function for processtrak interactive input of a frequency
# with read of default value with noprovenance. 
import sys
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv

   
if sys.version_info[0] < 3:
    input = raw_input   # backwards compatibility with Python 2.x
    pass



def enter_frequency(_xmldoc,_element,descr,defaulttagname):
    priorfreq=_xmldoc.xpathsinglecontext(_element,defaulttagname,default=None,noprovenance=True)
    if priorfreq is not None:
        priorfreq_val = numericunitsv.fromxml(_xmldoc,priorfreq).value('Hz')
        freq_text = input("Enter %s frequency in Hz (default %f): " % (descr,priorfreq_val))      
	if len(freq_text.strip())==0:
            freq = priorfreq_val
            pass
        else: 
	    freq=float(freq_text.strip())
	    pass
        pass
    else: 
        freq_text = input("Enter %s frequency in Hz: " % (descr))      
        freq=float(freq_text.strip())
	pass
    return freq
