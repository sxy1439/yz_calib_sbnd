import math 
import json
import ROOT
import numpy as np

# Open the ROOT file
root_file = ROOT.TFile("/exp/sbnd/app/users/kshiba/YZ_NonUniformity/srcs/yz_calib_sbnd/output_files/run_2_sample/run_2_gen_2_BNBLight_data_v10_14_02_02(200usEast80usWest).root", "READ")

# Get the TH2D histogram from the file
histogram_name = ["CzyHist_0_0",
                  "CzyHist_1_0",
                  "CzyHist_2_0",    
                  "CzyHist_0_1",
                  "CzyHist_1_1",
                  "CzyHist_2_1"]

pydict1 = {"AnodePlane:anode110": {}}
pydict2 = {"AnodePlane:anode120": {}}
pydict3 = {"AnodePlane:anode111": {}}
pydict4 = {"AnodePlane:anode121": {}}


for name in histogram_name:
    th2d_histogram = root_file.Get(name)

#    th2d_histogram = ROOT.TH2D(th2d_histogram)

    # Get the dimensions of the TH2D histogram
    n_bins_x = th2d_histogram.GetNbinsX()
    n_bins_y = th2d_histogram.GetNbinsY()
    
    # Create a NumPy array to store the histogram content
    histogram_array = np.array([])
    histogram_array = np.zeros((n_bins_x, n_bins_y))
    
    # Fill the NumPy array with the histogram content
    for i in range(1, n_bins_x + 1):
        for j in range(1, n_bins_y + 1):
            bin_content = th2d_histogram.GetBinContent(i, j)
            if math.isnan(bin_content): 
                bin_content = 1

            histogram_array[i - 1, j - 1] = str(round(bin_content, 3))

            # Close the ROOT file

            #  name_map("anode110")="East TPC, so anthing with x_0";
            #  name_map("anode120")="East TPC, so anthing with x_0";
            #  name_map("anode111")="West TPC, so anthing with x_1";
            #  name_map("anode121")="West TPC, so anthing with x_1";

    if name[-1] == "0":
        if "0_0" in name:
            pydict1["AnodePlane:anode110"]["0"] = histogram_array.tolist()
        if "1_0" in name:
            pydict1["AnodePlane:anode110"]["1"] = histogram_array.tolist()
        if "2_0" in name:
            pydict1["AnodePlane:anode110"]["2"] = histogram_array.tolist()
    if name[-1] == "0":
        if "0_0" in name:
            pydict2["AnodePlane:anode120"]["0"] = histogram_array.tolist()
        if "1_0" in name:
            pydict2["AnodePlane:anode120"]["1"] = histogram_array.tolist()
        if "2_0" in name:
            pydict2["AnodePlane:anode120"]["2"] = histogram_array.tolist()
    if name[-1] == "1":
        if "0_1" in name:
            pydict3["AnodePlane:anode111"]["0"] = histogram_array.tolist()
        if "1_1" in name:
            pydict3["AnodePlane:anode111"]["1"] = histogram_array.tolist()
        if "2_1" in name:
            pydict3["AnodePlane:anode111"]["2"] = histogram_array.tolist()
    if name[-1] == "1":
        if "0_1" in name:
            pydict4["AnodePlane:anode121"]["0"] = histogram_array.tolist()
        if "1_1" in name:
            pydict4["AnodePlane:anode121"]["1"] = histogram_array.tolist()
        if "2_1" in name:
            pydict4["AnodePlane:anode121"]["2"] = histogram_array.tolist()
    


# Merge all dictionaries into a single top-level dictionary
combined_dict = {}
for d in [pydict1, pydict2, pydict3, pydict4]:
    combined_dict.update(d)

# Write to a well-formed JSON file
with open("yzmap_gain_sbnd_run2.json", "w") as outfile:
    json.dump(combined_dict, outfile, indent=2)
 
root_file.Close()







"""
ORIGNAL CODE HERE. 


import math 
import json
import ROOT
import numpy as np

# Open the ROOT file
root_file = ROOT.TFile("run_2_gen_2_BNBLight_data_v10_14_02_02.root", "READ")

# Get the TH2D histogram from the file
histogram_name = {"TPCEE_P0_scale",
                  "TPCEE_P1_scale",
                  "TPCEE_P2_scale",
                  "TPCEW_P0_scale",
                  "TPCEW_P1_scale",
                  "TPCEW_P2_scale",
                  "TPCWE_P0_scale",
                  "TPCWE_P1_scale",
                  "TPCWE_P2_scale",
                  "TPCWW_P0_scale",
                  "TPCWW_P1_scale",
                  "TPCWW_P2_scale"}

pydict1 = {"AnodePlane:anode110": {}}
pydict2 = {"AnodePlane:anode120": {}}
pydict3 = {"AnodePlane:anode111": {}}
pydict4 = {"AnodePlane:anode121": {}}
pydict5 = {"AnodePlane:anode112": {}}
pydict6 = {"AnodePlane:anode122": {}}
pydict7 = {"AnodePlane:anode113": {}}
pydict8 = {"AnodePlane:anode123": {}}

for name in histogram_name:
    th2d_histogram = root_file.Get(name)

#    th2d_histogram = ROOT.TH2D(th2d_histogram)

    # Get the dimensions of the TH2D histogram
    n_bins_x = th2d_histogram.GetNbinsX()
    n_bins_y = th2d_histogram.GetNbinsY()
    
    # Create a NumPy array to store the histogram content
    histogram_array = np.array([])
    histogram_array = np.zeros((n_bins_x, n_bins_y))
    
    # Fill the NumPy array with the histogram content
    for i in range(1, n_bins_x + 1):
        for j in range(1, n_bins_y + 1):
            bin_content = th2d_histogram.GetBinContent(i, j)
            if math.isnan(bin_content): 
                bin_content = 1

            histogram_array[i - 1, j - 1] = str(round(bin_content, 3))

            # Close the ROOT file

            #  name_map("anode110")="EE";
            #  name_map("anode120")="EE";
            #  name_map("anode111")="EW";
            #  name_map("anode121")="EW";
            #  name_map("anode112")="WE";
            #  name_map("anode122")="WE";
            #  name_map("anode113")="WW";
            #  name_map("anode123")="WW";

    if "EE" in name:
        if "P0" in name:
            pydict1["AnodePlane:anode110"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict1["AnodePlane:anode110"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict1["AnodePlane:anode110"]["2"] = histogram_array.tolist()
    if "EE" in name:
        if "P0" in name:
            pydict2["AnodePlane:anode120"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict2["AnodePlane:anode120"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict2["AnodePlane:anode120"]["2"] = histogram_array.tolist()
    if "EW" in name:
        if "P0" in name:
            pydict3["AnodePlane:anode111"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict3["AnodePlane:anode111"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict3["AnodePlane:anode111"]["2"] = histogram_array.tolist()
    if "EW" in name:
        if "P0" in name:
            pydict4["AnodePlane:anode121"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict4["AnodePlane:anode121"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict4["AnodePlane:anode121"]["2"] = histogram_array.tolist()
    if "WE" in name:
        if "P0" in name:
            pydict5["AnodePlane:anode112"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict5["AnodePlane:anode112"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict5["AnodePlane:anode112"]["2"] = histogram_array.tolist()
    if "WE" in name:
        if "P0" in name:
            pydict6["AnodePlane:anode122"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict6["AnodePlane:anode122"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict6["AnodePlane:anode122"]["2"] = histogram_array.tolist()
    if "WW" in name:
        if "P0" in name:
            pydict7["AnodePlane:anode113"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict7["AnodePlane:anode113"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict7["AnodePlane:anode113"]["2"] = histogram_array.tolist()
    if "WW" in name:
        if "P0" in name:
            pydict8["AnodePlane:anode123"]["0"] = histogram_array.tolist()
        if "P1" in name:
            pydict8["AnodePlane:anode123"]["1"] = histogram_array.tolist()
        if "P2" in name:
            pydict8["AnodePlane:anode123"]["2"] = histogram_array.tolist()


# Merge all dictionaries into a single top-level dictionary
combined_dict = {}
for d in [pydict1, pydict2, pydict3, pydict4, pydict5, pydict6, pydict7, pydict8]:
    combined_dict.update(d)

# Write to a well-formed JSON file
with open("yzmap_gain_sbnd_run2.json", "w") as outfile:
    json.dump(combined_dict, outfile, indent=2)
 
root_file.Close()

"""