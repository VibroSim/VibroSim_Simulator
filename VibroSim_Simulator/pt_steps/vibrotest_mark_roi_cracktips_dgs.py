import sys
import subprocess
import ast
import numpy as np
import numbers

import numpy as np
from matplotlib import pyplot as pl

import dataguzzler as dg
from dataguzzler import metadata as dgm
import dg_file as dgf

from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.dc_value import stringvalue as stringv
from limatix.dc_value import hrefvalue as hrefv
from VibroSim_Simulator import format_modes

   
if sys.version_info[0] < 3:
    input = raw_input   # backwards compatibility with Python 2.x
    pass



def enter_multicoords(_xmldoc,_element):
    priorROIlowerleft=ast.literal_eval(_xmldoc.xpathsinglecontextstr(_element,"dc:dgs_roicorner1",default="None",noprovenance=True))
    priorROIupperright=ast.literal_eval(_xmldoc.xpathsinglecontextstr(_element,"dc:dgs_roicorner2",default="None",noprovenance=True))
    priorside1tip=ast.literal_eval(_xmldoc.xpathsinglecontextstr(_element,"dc:dgs_cracktipcoords_side1",default="None",noprovenance=True))
    priorside2tip=ast.literal_eval(_xmldoc.xpathsinglecontextstr(_element,"dc:dgs_cracktipcoords_side2",default="None",noprovenance=True))

    priors = ( priorROIlowerleft, priorROIupperright, priorside1tip, priorside2tip )
    
    default_string=", ".join([str(prior) for prior in priors])

    if default_string == "None, None, None, None":
        coord_text = input("Enter parenthesized pairs of coordinates,\nseparated by commas: ")      
        coords = ast.literal_eval(coord_text)
        
        pass
    else:
        coord_text = input("Enter parenthesized pairs of coordinates,\nseparated by commas (default %s): " % (default_string))      
        if len(coord_text.strip())==0:
            coords = priors
            pass
        else:
            coords = ast.literal_eval(coord_text)
            pass
        pass
    
    return coords  # ( ROIlowerleft, ROIupperright, side1tip, side2tip )

def run(_xmldoc,_element,dc_dgsfile_href,
        dc_exc_t0_numericunits,
        dc_exc_t3_numericunits,
        dc_simulationcameranetd_numericunits):  # Simulated camera NETD is used to generate same plot bounds as vibrosim_heatflow_analysis.py from crack_heatflow package
    
    print("DGS file: %s" % (dc_dgsfile_href.getpath()))
    
    print("Based on the displayed .dgs file or known consistent orientation,")
    print("enter coordinates of two ROI corners (lower left, upper right)")
    print("followed by side 1 (left or bottom) crack tip location, and ")
    print("crack side 2 (right or top) crack tip location.")
    print(" ")
    print("i.e. (roi_left,roi_bot),(roi_right,roi_top),(cracktip1_x,cracktip1_y),(cracktip2_x,cracktip2_y)")
    print(" ")

    print("If one side of the crack does not exist, click for that tip the")
    print("symmetric position about the crack center where that tip would be.")
    print(" ")
    print("Remember that the tip positions should represent the physical crack")
    print("tip locations, not the bounds of the heating. Use known crack size")
    print("if possible.")
    print(" ")
    print("Click the desired point in the scope, then middle mouse button")
    print("in this window to paste the coordinates")
    
    dgscope_proc = subprocess.Popen(["dg_scope_sa",dc_dgsfile_href.getpath()])

    # Test for good input
    good_input = False
    
    while not good_input:

        try:
            (ROIlowerleft,ROIupperright,side1tip,side2tip) = enter_multicoords(_xmldoc,_element)
            
            if len(ROIlowerleft) != 2 or not isinstance(ROIlowerleft[0],numbers.Number) or not isinstance(ROIlowerleft[1],numbers.Number):
                raise ValueError("Coordinates must be numbers")
            if len(ROIupperright) != 2 or not isinstance(ROIupperright[0],numbers.Number) or not isinstance(ROIupperright[1],numbers.Number):
                raise ValueError("Coordinates must be numbers")
            if len(side1tip) != 2 or not isinstance(side1tip[0],numbers.Number) or not isinstance(side1tip[1],numbers.Number):
                raise ValueError("Coordinates must be numbers")
            if len(side2tip) != 2 or not isinstance(side2tip[0],numbers.Number) or not isinstance(side2tip[1],numbers.Number):
                raise ValueError("Coordinates must be numbers")


            t0 = dc_exc_t0_numericunits.value("s") # expected welder start time                                                                                 
            t3 = dc_exc_t3_numericunits.value("s") # expected welder end time                                                                                   
            
            (dgsmetadata,wfmdict) = dgf.loadsnapshot(dc_dgsfile_href.getpath())
            DiffStack = wfmdict["DiffStack"]
            DiffStack_dx = dgm.GetMetaDatumWIDbl(DiffStack,"Step1",1.0)
            DiffStack_x = dgm.GetMetaDatumWIDbl(DiffStack,"IniVal1",0.0) + np.arange(DiffStack.data.shape[0],dtype='d')*DiffStack_dx
            DiffStack_dy = dgm.GetMetaDatumWIDbl(DiffStack,"Step2",1.0)
            DiffStack_y = dgm.GetMetaDatumWIDbl(DiffStack,"IniVal2",0.0) + np.arange(DiffStack.data.shape[1],dtype='d')*DiffStack_dy
            DiffStack_dt = dgm.GetMetaDatumWIDbl(DiffStack,"Step3",1.0)
            DiffStack_t = dgm.GetMetaDatumWIDbl(DiffStack,"IniVal3",0.0) + np.arange(DiffStack.data.shape[2],dtype='d')*DiffStack_dt
            DiffStack_units = dgm.GetMetaDatumWIStr(DiffStack,"AmplUnits","Volts")
            
            target_time = t0 + ((t3-t0)*(2.0/3.0)) # 2/3rds of way from t0 to t3                                                                                
            
            target_frameidx = np.argmin(abs(DiffStack_t-target_time))
            
            target_frame = DiffStack.data[:,:,target_frameidx]

            
            if DiffStack_units != "K":                                                                                                                          
                raise ValueError("DiffStack channel in %s has incorrect units (got \"%s\"; expected \"K\")." % (dc_dgsfile_href.getpath(),DiffStack_units))     
            
            camera_netd = dc_simulationcameranetd_numericunits.value("K")                                                                                       

            pl.figure()
            pl.imshow(target_frame.T,extent=(DiffStack_x[0]-DiffStack_dx/2.0,DiffStack_x[-1]+DiffStack_dx/2.0,DiffStack_y[0]-DiffStack_dy/2.0,DiffStack_y[-1]+DiffStack_dy/2.0),origin="lower",cmap="hot",vmin=-camera_netd,vmax=camera_netd*9)
            pl.plot(side1tip[0],side1tip[1],'x')
            pl.plot(side2tip[0],side2tip[1],'*')
            pl.plot((ROIlowerleft[0],ROIlowerleft[0],ROIupperright[0],ROIupperright[0],ROIlowerleft[0]),
                    (ROIlowerleft[1],ROIupperright[1],ROIupperright[1],ROIlowerleft[1],ROIlowerleft[1]),'-')
            pl.legend(("Side 1 (left/bottom) tip","Side 2 (right/top) tip","ROI"))
            pl.colorbar()
            pl.xlabel("%s (%s)" % (dgm.GetMetaDatumWIStr(DiffStack,"Coord1","X Position"),dgm.GetMetaDatumWIStr(DiffStack,"Units1","px")))
            pl.ylabel("%s (%s)" % (dgm.GetMetaDatumWIStr(DiffStack,"Coord2","Y Position"),dgm.GetMetaDatumWIStr(DiffStack,"Units2","px")))
            pl.title("t = %f s" % (DiffStack_t[target_frameidx]))

            print("Now view and assess the extracted thermal image with locations marked")
            print("Then close the thermal image window.")
            pl.show()

            OK = input("Is this OK [y/N]: ")     
            if OK.strip().lower()=="y" or OK.strip().upper()=="yes":
                good_input = True
                pass
            else:
                print("Retrying...")
                pass
            pass
        except SyntaxError as e:
            print("Syntax error: %s; try again" % (str(e)))
            pass
        except ValueError as e:
            print("Value error: %s; try again" % (str(e)))
            pass
        pass
    pass

    print("Now close dg_scope if you haven't already.")

    dgscope_proc.communicate() # Wait for COMSOL to finish

    return { 
        "dc:dgs_roicorner1": stringv(ROIlowerleft), 
        "dc:dgs_roicorner2": stringv(ROIupperright), 
        "dc:dgs_cracktipcoords_side1": stringv(side1tip), 
        "dc:dgs_cracktipcoords_side2": stringv(side2tip), 
    }
