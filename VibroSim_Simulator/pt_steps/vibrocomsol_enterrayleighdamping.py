import sys
import subprocess
import traceback
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_generic # enter_frequency
from VibroSim_Simulator import open_via_helper

def run(_xmldoc,_element):

    # format_modes.print_comsol_modes(dc_modalfreqs_href.getpath())

    # print("\n")

    # print("Based on the modal analysis in COMSOL (under")
    # print("results/vibro_modal_plot)")

    # imagehelper_wait=lambda: None
    # image_viewer_text=None

    # if dc_sweep_spectrum_href is not None: 
    #     print("and the spectrum plot")
    #     image_viewer_text = "and the image viewer "
    #     imagehelper_wait = open_via_helper.open_via_helper(dc_sweep_spectrum_href.getpath())
    #     pass

    # print("identify a frequency for the burst. ")
    print('Enter damping factors.')

    # comsolproc = subprocess.Popen(["comsol","-open",dc_modalcalc_comsol_href.getpath()])

    # Extract prior value of dc:excitation_frequency with noprovenance as default
    # Check that a number was given
    number_input=False
    while not number_input:
        try:
            alpha = enter_generic.enter_generic(_xmldoc,_element,"default","dc:spcrayleighdamping_alpha","s^-1")
            number_input=True
        except ValueError:
            pass

    number_input=False
    while not number_input:
        try:
            beta= enter_generic.enter_generic(_xmldoc,_element,"default","dc:spcrayleighdamping_beta","s")
            number_input=True
        except ValueError:
            pass

    # print("Now close COMSOL %sif you haven't already." % (image_viewer_text))
    # comsolproc.communicate()
    # imagehelper_wait()

    return { "dc:spcrayleighdamping_alpha": numericunitsv(alpha,"s^-1") , "dc:spcrayleighdamping_beta": numericunitsv(beta,"s") }
