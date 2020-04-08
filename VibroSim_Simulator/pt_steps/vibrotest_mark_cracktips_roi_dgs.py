import sys
import subprocess
import ast
import numpy as np
import numbers

from limatix.dc_value import numericunitsvalue as numericunitsv
from limatix.dc_value import stringvalue as stringv
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
    
    default_string=", ".join(priors)

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

def run(_xmldoc,_element,dc_dgsfile_href):

    print("DGS file: %s" % (dc_dgsfile_href.getpath()))
    
    print("Based on the displayed .dgs file or known consistent orientation,")
    print("enter coordinates of two ROI corners (lower left, upper right)")
    print("followed by side 1 (left or bottom) crack tip location, and ")
    print("crack side 2 (right or top) crack tip location.")
    print(" ")

    print("If one side of the crack does not exist, click for that tip the")
    print("symmetric position about the crack center where that tip would be.")
    print(" ")
    print("Click the desired point in the scope, then middle mouse button")
    print("in this window to paste the coordinates")
    
    dgscope_proc = subprocess.Popen(["dc_scope_sa",dc_dgsfile_href.getpath()])

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
                
            
            good_input=True
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
