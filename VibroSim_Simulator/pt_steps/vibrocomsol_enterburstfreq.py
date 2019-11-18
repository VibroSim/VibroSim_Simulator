import sys
import traceback
import numpy as np

from limatix.dc_value import numericunitsvalue as numericunitsv

from VibroSim_Simulator import format_modes
from VibroSim_Simulator import enter_frequency

def run(_xmldoc,_element, dc_modalfreqs_href):

    format_modes.print_comsol_modes(dc_modalfreqs_href.getpath())

    print("\n")

    # Extract prior value of dc:excitation_frequency with noprovenance as default
    freq = enter_frequency.enter_frequency(_xmldoc,_element,"burst","dc:excitation_frequency")

    return { "dc:excitation_frequency": numericunitsv(freq,"Hz") }
