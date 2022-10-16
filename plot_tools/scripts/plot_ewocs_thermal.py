import sys
import numpy as np

from ewoc_cmdln import ewoc_parser
from ewoc_utils import plot_EWOC_by_frag_temps

# Get arguments
kwargs=ewoc_parser.parse_args()

plot_EWOC_by_frag_temps(load=True,
                      print_every_n=500,
                      # hist_weight=lambda x: x*(1-x),
                      **vars(kwargs))

