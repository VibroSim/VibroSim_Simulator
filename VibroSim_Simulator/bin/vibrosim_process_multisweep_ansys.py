import os
import os.path
import sys

import numpy as np
from matplotlib import pyplot as pl

try:
    # py2.x
    from urllib import pathname2url
    from urllib import url2pathname
    from urllib import quote
    from urllib import unquote
    pass
except ImportError:
    # py3.x
    from urllib.request import pathname2url
    from urllib.request import url2pathname
    from urllib.parse import quote
    from urllib.parse import unquote
    pass


from limatix.dc_value import hrefvalue as hrefv
from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator.process_multisweep import process_multisweep_ansys


def main(args=None):

    contactpoint_name = "Horn Contact Point"
    
    
    if args is None:
        args=sys.argv
        pass

    if len(args) < 8:
        print("usage: vibrosim_process_multisweep_ansys <laser_name> <crack_name> <plotdir> <endcrop_seconds> <seg1_table.xlsx> <seg2_table.xlsx> <seg3_table.xlsx> <seg4_table.xlsx> <output_filename.csv.bz2> ")
        print(" ")
        print("Example: vibrosim_process_multisweep_ansys \"SS TE Laser Spot\" \"SN10 SS Tip Crack\" .  0.05 \"Frequency Response 0Hz to 9990 Hz.xlsx\" \"Frequency Response 10000Hz to 99900Hz.xlsx\" \"Frequency Response 100000Hz to 999000Hz.xlsx\" specimen_response.csv.bz2")
        sys.exit(0)
        pass

    laser_name=args[1]
    crack_name=args[2]

    plotdir = args[3]
    endcrop = float(args[4])

    
    #seg1_table = args[5]
    #seg2_table = args[6]
    #seg3_table = args[7]
    #seg4_table = args[8]

    seg_tables = args[5:-1]

    output_filename=args[-1]

    if not(output_filename.endswith(".csv.bz2")) and not(output_filename.endswith(".csv")):
        print("Output filename \"%s\" should have a .csv or .csv.bz2 extension" % (output_filename))
        sys.exit(1)
        pass
    
    
    dt=1e-6
    impulseexcitation_width=2.0e-6

    (output_filename_ref, plotpaths) = process_multisweep_ansys(laser_name,crack_name,plotdir,"",dt,impulseexcitation_width,endcrop,seg_tables, output_filename)
    
    
    pass
