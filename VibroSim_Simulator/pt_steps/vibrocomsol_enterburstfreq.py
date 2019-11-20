import sys
import subprocess
import traceback
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_frequency
from VibroSim_Simulator import open_via_helper

def run(_xmldoc,_element, dc_modalfreqs_href,dc_modalcalc_comsol_href,dc_sweep_spectrum_href=None):

    format_modes.print_comsol_modes(dc_modalfreqs_href.getpath())

    print("\n")

    print("Based on the modal analysis in COMSOL (under")
    print("results/vibro_modal_plot)")

    imagehelper_wait=lambda: None
    image_viewer_text=None

    if dc_sweep_spectrum_href is not None: 
        print("and the spectrum plot")
        image_viewer_text = "and the image viewer "
        imagehelper_wait = open_via_helper.open_via_helper(dc_sweep_spectrum_href.getpath())
        pass

    print("identify a frequency for the burst. ")

    comsolproc = subprocess.Popen(["comsol","-open",dc_modalcalc_comsol_href.getpath()])

    # Extract prior value of dc:excitation_frequency with noprovenance as default
    freq = enter_frequency.enter_frequency(_xmldoc,_element,"burst","dc:excitation_frequency")


    print("Now close COMSOL %sif you haven't already." % (image_viewer_text))
    comsolproc.communicate()
    imagehelper_wait()


    return { "dc:excitation_frequency": numericunitsv(freq,"Hz") }
