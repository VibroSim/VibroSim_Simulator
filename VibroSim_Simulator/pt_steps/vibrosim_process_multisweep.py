import os
import os.path
import sys

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

from VibroSim_Simulator.process_multisweep import process_multisweep_from_files
from VibroSim_Simulator.function_as_script import scriptify


def run(_xmldoc,_element,
        dc_dest_href,
        dc_measident_str,
        dc_seg1_xducercontactprobe_displ_href,
        dc_seg1_xducercontactprobe_vel_href,
        dc_seg1_laser_displ_href,
        dc_seg1_laser_vel_href,
        dc_seg1_crackcenterstress_href,
        dc_seg2_xducercontactprobe_displ_href,
        dc_seg2_xducercontactprobe_vel_href,
        dc_seg2_laser_displ_href,
        dc_seg2_laser_vel_href,
        dc_seg2_crackcenterstress_href,
        dc_seg3_xducercontactprobe_displ_href,
        dc_seg3_xducercontactprobe_vel_href,
        dc_seg3_laser_displ_href,
        dc_seg3_laser_vel_href,
        dc_seg3_crackcenterstress_href,
        dc_seg4_xducercontactprobe_displ_href,
        dc_seg4_xducercontactprobe_vel_href,
        dc_seg4_laser_displ_href,
        dc_seg4_laser_vel_href,
        dc_seg4_crackcenterstress_href,
        dc_endcrop_numericunits=numericunitsv(250e-3,"s")
):

    endcrop = dc_endcrop_numericunits.value('s')
    dt=1e-6
    impulseexcitation_width=2.0e-6

    xducer_velspec_filepaths = [ 
        dc_seg1_xducercontactprobe_vel_href.getpath(),
        dc_seg2_xducercontactprobe_vel_href.getpath(),
        dc_seg3_xducercontactprobe_vel_href.getpath(),
        dc_seg4_xducercontactprobe_vel_href.getpath(),
    ]

    xducer_displspec_filepaths = [ 
        dc_seg1_xducercontactprobe_displ_href.getpath(),
        dc_seg2_xducercontactprobe_displ_href.getpath(),
        dc_seg3_xducercontactprobe_displ_href.getpath(),
        dc_seg4_xducercontactprobe_displ_href.getpath(),
    ]

    laser_velspec_filepaths = [ 
        dc_seg1_laser_vel_href.getpath(),
        dc_seg2_laser_vel_href.getpath(),
        dc_seg3_laser_vel_href.getpath(),
        dc_seg4_laser_vel_href.getpath(),
    ]

    laser_displspec_filepaths = [ 
        dc_seg1_laser_displ_href.getpath(),
        dc_seg2_laser_displ_href.getpath(),
        dc_seg3_laser_displ_href.getpath(),
        dc_seg4_laser_displ_href.getpath(),
    ]

    crackcenterstressspec_filepaths = [ 
        dc_seg1_crackcenterstress_href.getpath(),
        dc_seg2_crackcenterstress_href.getpath(),
        dc_seg3_crackcenterstress_href.getpath(),
        dc_seg4_crackcenterstress_href.getpath(),
    ]
    
    output_href=hrefv(quote(dc_measident_str+"_dynamicmodel.csv.bz2"),dc_dest_href)


    # Can use scriptify(process_multisweep_from_files) in next line to debug processing script
    (output_filename, plotpaths) = process_multisweep_from_files(xducer_velspec_filepaths,
                                                                 xducer_displspec_filepaths,
                                                                 laser_velspec_filepaths,
                                                                 laser_displspec_filepaths,
                                                                 crackcenterstressspec_filepaths,
                                                                 output_href.getpath(),
                                                                 plotdir=dc_dest_href.getpath(),
                                                                 plotprefix=dc_measident_str,
                                                                 num_segments=4,
                                                                 dt=dt,
                                                                 impulseexcitation_width=impulseexcitation_width, 
                                                                 endcrop=endcrop)
    ret = {
        "dc:dynamicmodel": output_href,
    }

    for key in plotpaths:
        ret["dc:%s" % (key)] = hrefv(pathname2url(plotpaths[key]),".")
        pass
    return ret
    
