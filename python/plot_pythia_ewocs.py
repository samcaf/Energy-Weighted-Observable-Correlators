import sys
import numpy as np
from ewoc_utils import *

plot_EWOC_by_sub_rads(sys.argv)

# return
# # Get command line input
# process_str = sys.argv[1]

# file_prefix = ewoc_file_label(sys.argv)

# data_dir = file_prefix + '.txt'
# hist_dir = file_prefix + '.npz'

# try:
#     # Loading data if it exists
#     loaded_dict = np.load(hist_dir, allow_pickle=True)
#     hist_dict = {}
#     for key in loaded_dict.keys():
#         hist_dict[key] = loaded_dict[key][()]
#     del loaded_dict

# except FileNotFoundError as e:
#     # Otherwise, creating the data
#     print(e.message())
#     data_dict = ewoc_text_to_dict(data_dir)
#     hist_dict = ewoc_dict_to_hists(data_dict) 
#     np.savez(hist_dir, **hist_dict)

# plot_hist_dict(hist_dict)
