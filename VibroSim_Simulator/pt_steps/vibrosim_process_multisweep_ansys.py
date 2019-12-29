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
        dc_seg1_table_href,
        dc_seg2_table_href,
        dc_seg3_table_href,
        dc_seg4_table_href,        
        dc_endcrop_numericunits=numericunitsv(250e-3,"s")
):

    endcrop = dc_endcrop_numericunits.value('s')
    dt=1e-6
    impulseexcitation_width=2.0e-6
    
    
    output_href=hrefv(quote(dc_measident_str+"_dynamicmodel.csv.bz2"),dc_dest_href)

    
    # Can use scriptify(process_multisweep_ansys) in next line to debug processing script
    (output_filename, plotpaths) = process_multisweep_ansys(dc_laser_name_str,
                                                            dc_crack_name_str,
                                                            dc_dest_href.getpath(),
                                                            dt,
                                                            impulseexcitation_width,
                                                            endcrop,
                                                            [ dc_seg1_table_href.getpath(),
                                                              dc_seg2_table_href.getpath(),
                                                              dc_seg3_table_href.getpath(),
                                                              dc_seg4_table_href.getpath() ],
                                                            output_href.getpath())
    ret = {
        "dc:dynamicmodel": output_href,
    }

    for key in plotpaths:
        ret["dc:%s" % (key)] = hrefv(pathname2url(plotpaths[key]),".")
        pass
    return ret
    
