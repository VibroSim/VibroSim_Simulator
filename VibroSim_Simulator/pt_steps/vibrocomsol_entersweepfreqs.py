import sys
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv
from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_frequency


def run(_xmldoc,_element, dc_modalfreqs_href):

    format_modes.print_comsol_modes(dc_modalfreqs_href.getpath())


    print("\n")


    startfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep start","dc:sweep_start_frequency")
    stepfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep step","dc:sweep_step_frequency")
    endfreq = enter_frequency.enter_frequency(_xmldoc,_element,"sweep end","dc:sweep_end_frequency")

    return { 
      "dc:sweep_start_frequency": numericunitsv(startfreq,"Hz"), 
      "dc:sweep_step_frequency": numericunitsv(stepfreq,"Hz"), 
      "dc:sweep_end_frequency": numericunitsv(endfreq,"Hz"), 
    }
