import sys
import subprocess
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv
from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_frequency


def run(_xmldoc,_element, dc_modalfreqs_href,dc_modalcalc_comsol_href):

    format_modes.print_comsol_modes(dc_modalfreqs_href.getpath())


    print("\n")

    print("Based on the modal analysis in COMSOL (under")
    print("results/vibro_modal_plot), identify a band of frequencies")
    print("(start,step,end) for the sweep. ")

    comsolproc = subprocess.Popen(["comsol","-open",dc_modalcalc_comsol_href.getpath()])

    startfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep start","dc:sweep_start_frequency")
    stepfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep step","dc:sweep_step_frequency")
    endfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep end","dc:sweep_end_frequency")

    print("Now close COMSOL if you haven't already.")

    comsolproc.communicate() # Wait for COMSOL to finish

    return { 
      "dc:sweep_start_frequency": numericunitsv(startfreq,"Hz"), 
      "dc:sweep_step_frequency": numericunitsv(stepfreq,"Hz"), 
      "dc:sweep_end_frequency": numericunitsv(endfreq,"Hz"), 
    }
