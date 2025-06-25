#include "file_handling.h"
#include <iostream>
#include <cstdlib>  // for std::exit

std::array<float, 2> getELifetime(const std::string& cintyp){
  if (cintyp == "doMC2023B_sub1" || cintyp == "doMC2023B_full") {
    std::array<float, 2> e_callife = {10.e3, 10.e3};// in us
    return e_callife;
  } else if(cintyp == "doData_dev" || cintyp == "doData_dev_scaling1" || cintyp == "doData_dev_ratio1"){
    std::array<float, 2> e_callife = {60.52e3, 92.46e3};// in us
    return e_callife;
  }
  else {
    std::array<float, 2> e_callife = {100.e3, 100.e3};// in us
    return e_callife;
  }
}

TString getOutputNameChannelPerTheta(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "output_files/channel_per_thetaxz_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/channel_per_thetaxz_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_per_thetaxz_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_per_thetaxz_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/channel_per_thetaxz_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/channel_per_thetaxz_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/channel_per_thetaxz_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/channel_per_thetaxz_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/channel_per_thetaxz_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameChannelByAngle(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/channel_by_angle_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/channel_by_angle_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_by_angle_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_by_angle_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/channel_by_angle_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/channel_by_angle_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/channel_by_angle_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/channel_by_angle_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/channel_by_angle_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/channel_by_angle_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameChannel(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/channel_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/channel_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/channel_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/channel_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/channel_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/channel_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/channel_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/channel_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameAng(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/ang_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/ang_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/ang_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/ang_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/ang_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/ang_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/ang_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/ang_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/ang_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/ang_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameYZ(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/yz_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/yz_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/yz_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/yz_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/yz_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/yz_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/yz_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/yz_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/yz_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/yz_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameX(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/x_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/x_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/x_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/x_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/x_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/x_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/x_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/x_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/x_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/x_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNamedQdx(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/dQdx_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/dQdx_data_dev.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/dQdx_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/dQdx_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/dQdx_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/dQdx_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/dQdx_mcp2025Afeb.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/dQdx_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/dQdx_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/dQdx_mc2024B_full.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameOfflinePitch(const std::string& cintyp) {
  if (cintyp == "doData_dev") {
    return "../output_files/offline_pitch_cal_data_dev.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/offline_pitch_cal_mcp2025av3.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}




std::ifstream openInputFile(const std::string& cintyp) {
  std::ifstream fl;
  if (cintyp == "doData_runs17742_to_87_sub1") {
    fl.open("../input_files/data_runs17742_to_87_sub1.txt");
  } else if (cintyp == "doData_dev") {
    fl.open("../input_files/data2025av3_dev_reco2.txt");
  } else if (cintyp == "doData_dev_scaling1") {
    //fl.open("../input_files/dev_data_scaling1_lynn.txt");
    fl.open("../input_files/dev_data_scaling1th25_lynn.txt");
  } else if (cintyp == "doData_dev_ratio1") {
    //fl.open("../input_files/dev_data_ratio1_lynn.txt");
    fl.open("../input_files/dev_data_datath25_lynn.txt");
  } else if (cintyp == "doMC2025Av3") {
    fl.open("../input_files/mcp2025av3.txt");
  } else if (cintyp == "doMC2025Adebug") {
    fl.open("/exp/sbnd/data/users/lynnt/wirecell/debug_v10/hists_reco2.list");
  } else if (cintyp == "doMC2025A_Feb") {
    fl.open("../input_files/mcp2025a_febv10_04_03.txt");
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    fl.open("../input_files/reco2_data_runs/com_runs17742_to_87_ctcMap_sub1.txt");
  } else if (cintyp == "doMC2024B_sub123") {
    //fl.open("../sbnd_MCP2024B_ntuples/file_sbnd_MCP2024B_sub123.txt");
    fl.open("../input_files/sbnd_MCP2024B_ntuples/file_sbnd_MCP2024B_sub1.txt");
  } else if (cintyp == "doMC2024B_full") {
    fl.open("../input_files/sbnd_MCP2024B_ntuples/file_sbnd_MCP2024B.txt");
  } else {
    std::cerr << "Failed to open the file for " << cintyp << "\n";
    exit(1);
  }
  return fl;
}


