# enter_generic.py: Convenience function for processtrak interactive input
# with read of default value with noprovenance. 
import sys
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv

   
if sys.version_info[0] < 3:
    input = raw_input   # backwards compatibility with Python 2.x
    pass



def enter_generic(_xmldoc,_element,descr,defaulttagname,unit,label='Input:'):
    prior=_xmldoc.xpathsinglecontext(_element,defaulttagname,default=None,noprovenance=True)
    if prior is not None:
        prior_val = numericunitsv.fromxml(_xmldoc,prior).value(unit)
        val_text = input("Enter {} in {} (default {}): ".format(descr,unit,prior_val))      
        if len(val_text.strip())==0:
            val = prior_val
            pass
        else: 
            val=float(val_text.strip())
            pass
        pass
    else: 
        val_text = input("Enter {} in {}: ".format(descr,unit))      
        val=float(val_text.strip())
        pass
    return val
