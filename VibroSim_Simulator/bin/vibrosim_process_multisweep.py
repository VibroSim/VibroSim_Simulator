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

from VibroSim_Simulator.process_multisweep import process_multisweep


def main(args=None):

    if args is None:
        args=sys.argv
        pass

    if len(args) < 8:
        print("usage: vibrosim_process_multisweep <datadir_path>  <xducer_velspec_filename_template> <xducer_displspec_filename_template> <laser_velspec_filename_template> <laser_displspec_filename_template> <crackcenterstrainspec_filename_template> [output_filename.py.bz2] [plotdir] [endcrop_seconds]")
        print(" ")
        print("Example: vibrosim_process_multisweep .  meas1_withsegboundaries_xducercontactprobe_velspec_seg%d.txt meas1_withsegboundaries_xducercontactprobe_displspec_seg%d.txt meas1_withsegboundaries_laser_velspec_seg%d.txt meas1_withsegboundaries_laser_displspec_seg%d.txt meas1_withsegboundaries_crackcenterstrainspec_seg%d.txt specimen_modeling_out.py.bz2 . 0.1")
        print(" ")
        print("Optional templates (such as laser displacement, laser velocity, and crack center strain) can be specified as \"None\".")
        sys.exit(0)
        pass

    path=args[1]
    xducer_velspec_filename_template=args[2]
    xducer_displspec_filename_template=args[3]

    laser_velspec_filename_template=args[4]
    if laser_velspec_filename_template=="None":
        laser_velspec_filename_template = None
        pass

    laser_displspec_filename_template=args[5]
    if laser_displspec_filename_template=="None":
        laser_displspec_filename_template = None
        pass

    crackcenterstrainspec_filename_template=args[6]
    if crackcenterstrainspec_filename_template=="None":
        crackcenterstrainspec_filename_template = None
        pass

    if len(args) > 7:
        output_filename=args[7]
        pass
    else:
        output_filename=os.path.join(tempfile.tempdir(),"specimen_modeling_out.py.bz2")
        pass

    if len(args) > 8:
        plotdir = args[8]
        pass

    if len(args) > 9:
        endcrop = float(args[9])
        pass

    
    dt=1e-6
    impulseexcitation_width=2.0e-6

    

    (output_filename_ref, plotpaths) = process_multisweep(path,
                                                      xducer_velspec_filename_template,
                                                      xducer_displspec_filename_template,
                                                      laser_velspec_filename_template,
                                                      laser_displspec_filename_template,
                                                      crackcenterstrainspec_filename_template,
                                                      output_filename,
                                                      plotdir=plotdir,
                                                      plotprefix="",
                                                      num_segments=4,
                                                      dt=1e-6,
                                                      impulseexcitation_width=impulseexcitation_width, # (s) COMSOL parameter -- width of Gaussian pulse
                                                      endcrop=endcrop)
    pass
