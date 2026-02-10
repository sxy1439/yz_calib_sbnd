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
    return "../output_files/channel_per_thetaxz_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/channel_per_thetaxz_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/channel_per_thetaxz_data1e20.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_per_thetaxz_data_dev_scaling1.root";
  } else if (cintyp == "doDatadebug_u2v2w25deconlim") {
    return "../output_files/channel_per_thetaxz_data_debug_u2v2w25deconlim.root";
  }  else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_per_thetaxz_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/channel_per_thetaxz_mcp2025b5e18.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/channel_per_thetaxz_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025Adebug_lowth") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug_lowth.root";
  } else if (cintyp == "doMC2025Adebug_u35v3w25") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug_u35v3w25.root";
  } else if (cintyp == "doMC2025Adebug_u2v2w25deconlim") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug_u2v2w25deconlim.root";
  } else if (cintyp == "doMC2025Adebug_pgun") {
    return "../output_files/channel_per_thetaxz_mcp2025Adebug_pgun.root";
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
  } else if (cintyp == "doData1e20") {
    return "../output_files/channel_by_angle_data1e20.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_by_angle_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_by_angle_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/channel_by_angle_mcp2025b5e18.root";
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
  } else if (cintyp == "doData_fall25valII") {
    return "../output_files/calib_paper_output/channel_data_fall25valII.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/channel_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/channel_data1e20.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/channel_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/channel_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/channel_mcp2025b5e18.root";
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
  } else if (cintyp == "doMC2025C_fallvalII") {
    return "../output_files/calib_paper_output/channel_mc2025C_fallvalII.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameAng(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/ang_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_fall25valII") {
    return "../output_files/calib_paper_output/ang_data_fall25valII.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/ang_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/ang_data1e20.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/ang_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/ang_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/ang_mcp2025b5e18.root";
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
  } else if (cintyp == "doMC2025C_fallvalII") {
    return "../output_files/calib_paper_output/ang_mc2025C_fallvalII.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameYZ(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/yz_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_fall25valII") {
    return "../output_files/calib_paper_output/yz_data_fall25valII.root";
  } else if (cintyp == "doData_v10_14_02") {
    return "../output_files/2026Janfallprod_maps/yz_data2025_v10_14_02.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/yz_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/yz_data1e20.root";
  } else if (cintyp == "doData1e20_validation") {
    return "../output_files/yz_data1e20_validation.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/yz_data_dev_scaling1.root";
  } else if (cintyp == "doDatadebug_u2v2w25deconlim") {
    return "../output_files/yz_data_debug_u2v2w25deconlim.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/yz_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/yz_mcp2025b5e18.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/yz_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/yz_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025Adebug_lowth") {
    return "../output_files/yz_mcp2025Adebug_lowth.root";
  } else if (cintyp == "doMC2025Adebug_u35v3w25") {
    return "../output_files/yz_mcp2025Adebug_u35v3w25.root";
  } else if (cintyp == "doMC2025Adebug_u2v2w25deconlim") {
    return "../output_files/yz_mcp2025Adebug_u2v2w25deconlim.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/yz_mcp2025Afeb.root";    
  } else if (cintyp == "doMC2025B_Oct") {
    return "../output_files/yz_mcp2025b_OctFall_workshop.root";
  } else if (cintyp == "doData_dev25B_Oct") {
    return "../output_files/yz_data_dev25b_OctFall_workshop.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/yz_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/yz_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/yz_mc2024B_full.root";
  } else if (cintyp == "doMC2025C_fallvalII") {
    return "../output_files/calib_paper_output/yz_mc2025C_fallvalII.root";
  } else if (cintyp == "doMC2025_v10_14_02") {
    return "../output_files/2026Janfallprod_maps/yz_mc2025_v10_14_02.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNameX(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/x_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_fall25valII") {
    return "../output_files/calib_paper_output/x_data_fall25valII.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/x_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/x_data1e20.root";
  } else if (cintyp == "doData1e20_validation") {
    return "../output_files/x_data1e20_validation.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/x_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/x_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/x_mcp2025b5e18.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/x_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/x_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/x_mcp2025Afeb.root";
  } else if (cintyp == "doMC2025B_Oct") {
    return "../output_files/x_mcp2025b_OctFall_workshop.root";
  } else if (cintyp == "doData_dev25B_Oct") {
    return "../output_files/x_data_dev25b_OctFall_workshop.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/x_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/x_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/x_mc2024B_full.root";
  } else if (cintyp == "doMC2025C_fallvalII") {
    return "../output_files/calib_paper_output/x_mc2025C_fallvalII.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNamedQdx(const std::string& cintyp) {
  if (cintyp == "doData_runs17742_to_87_sub1") {
    return "../output_files/dQdx_data_runs17742_to_87_sub1.root";
  } else if (cintyp == "doData_fall25valII") {
    return "../output_files/calib_paper_output/dQdx_data_fall25valII.root";
  } else if (cintyp == "doData_dev") {
    return "../output_files/dQdx_data_dev.root";
  } else if (cintyp == "doData1e20") {
    return "../output_files/dQdx_data1e20.root";
  } else if (cintyp == "doData1e20_validation") {
    return "../output_files/dQdx_data1e20_validation.root";
  } else if (cintyp == "doData_dev_scaling1") {
    return "../output_files/dQdx_data_dev_scaling1.root";
  } else if (cintyp == "doData_dev_ratio1") {
    return "../output_files/dQdx_data_dev_ratio1.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/dQdx_mcp2025b5e18.root";
  } else if (cintyp == "doMC2025Av3") {
    return "../output_files/dQdx_mcp2025av3.root";
  } else if (cintyp == "doMC2025Adebug") {
    return "../output_files/dQdx_mcp2025Adebug.root";
  } else if (cintyp == "doMC2025A_Feb") {
    return "../output_files/dQdx_mcp2025Afeb.root";
  } else if (cintyp == "doMC2025B_Oct") {
    return "../output_files/dQdx_mcp2025b_OctFall_workshop.root";
  } else if (cintyp == "doData_dev25B_Oct") {
    return "../output_files/dQdx_data_dev25b_OctFall_workshop.root";
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    return "../output_files/dQdx_data_runs17742_to_87_ctcMap_sub1.root";
  } else if (cintyp == "doMC2024B_sub123") {
    return "../output_files/dQdx_mc2024B_sub123.root";
  } else if (cintyp == "doMC2024B_full") {
    return "../output_files/dQdx_mc2024B_full.root";
  } else if (cintyp == "doMC2025C_fallvalII") {
    return "../output_files/calib_paper_output/dQdx_mc2025C_fallvalII.root";
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
  } else if (cintyp == "doData_fall25valII") {
    fl.open("../input_files/calib_paper_samples/files_data25C_FallValidationII.txt");
  } else if (cintyp == "doData_v10_14_02") {
    fl.open("../input_files/2026Janfallprod_maps/files_data25_v10_14_02.txt");
  } else if (cintyp == "doData_dev") {
    fl.open("../input_files/data2025av3_dev_reco2.txt");
  } else if (cintyp == "doData1e20") {    // used for the yz maps produced to correct for in sbndcode
    fl.open("../input_files/data_1e20_sub11.txt");
  } else if (cintyp == "doData1e20_validation") {
    fl.open("../input_files/data_1e20_sub1_validation.txt");
  } else if (cintyp == "doData_dev_scaling1") {
    fl.open("../input_files/dev_data_uvwdeconlim.txt");
  } else if (cintyp == "doDatadebug_u2v2w25deconlim") {
    fl.open("../input_files/data_u2v2w25deconlim.txt");
  } else if (cintyp == "doData_dev_ratio1") {
    fl.open("../input_files/dev_data_datath25_lynn.txt");
  } else if (cintyp == "doMC2025B5e18") {
    fl.open("../input_files/mcp2025b5e18.txt");
  } else if (cintyp == "doMC2025Av3") {
    fl.open("../input_files/mcp2025av3.txt");
  } else if (cintyp == "doMC2025Adebug") {
    fl.open("../input_files/mcp2025av3_debug_jul24.txt");
  } else if (cintyp == "doMC2025Adebug_lowth") {
    fl.open("../input_files/mcp2025av3_debug_lowth.txt");
  } else if (cintyp == "doMC2025Adebug_u35v3w25") {
    fl.open("../input_files/mcp2025av3_debug_u35v3w25.txt");
  } else if (cintyp == "doMC2025Adebug_u2v2w25deconlim") {
    fl.open("../input_files/mc_u2v2w25deconlim.txt");
  } else if (cintyp == "doMC2025Adebug_pgun") {
    fl.open("../input_files/mcp2025av3_debug_pgun_nocathodecoss.txt");
  } else if (cintyp == "doMC2025A_Feb") {
    fl.open("../input_files/mcp2025a_febv10_04_03.txt");
  } else if (cintyp == "doMC2025B_Oct") {
    fl.open("../input_files/mcp2025b_fall_workshop.txt");
  } else if (cintyp == "doData_dev25B_Oct") {
    fl.open("../input_files/dev_data25b_fall_workshop.txt");
  } else if (cintyp == "doData_runs17742_to_87_ctcMap_sub1") {
    fl.open("../../reco2_data_runs/com_runs17742_to_87_ctcMap_sub1.txt");
  } else if (cintyp == "doMC2024B_sub123") {
    fl.open("../../sbnd_MCP2024B_ntuples/file_sbnd_MCP2024B_sub1.txt");
  } else if (cintyp == "doMC2024B_full") {
    fl.open("../../sbnd_MCP2024B_ntuples/file_sbnd_MCP2024B.txt");
  } else if (cintyp == "doMC2025C_fallvalII") {
    fl.open("../input_files/calib_paper_samples/files_MCP2025C_FallValidationII_subOneTenth.txt");
  } else if (cintyp == "doMC2025_v10_14_02") {
    fl.open("../input_files/2026Janfallprod_maps/files_MCP2025_v10_14_02_sub1.txt");
  } else {
    std::cerr << "Failed to open the file for " << cintyp << "\n";
    exit(1);
  }
  return fl;
}


std::ifstream openInputFile_official(const std::string& cintyp) {
  std::ifstream fl;
  if (cintyp == "doData1e20") {
    //fl.open("../input_files/data_1e20_sub1.txt");
    fl.open("../input_files/data_1e20_onefile.txt");
  } else if (cintyp == "doMC2025B5e18") {
    fl.open("../input_files/mcp2025b5e18.txt");
    std::cout<<"filse in use are from ../input_files/mcp2025b5e18.txt"<<std::endl;
  } else {
    std::cerr << "Failed to open the file for " << cintyp << "\n";
    exit(1);
  }
  return fl;
}

std::ifstream openInputFile_validation(const std::string& cintyp) {
  std::ifstream fl;
  if (cintyp == "doData1e20") {
    //fl.open("../input_files/data_1e20_sub1_validation.txt");
    fl.open("../input_files/data_1e20_onefile_validation.txt");
  } else if (cintyp == "doMC2025B5e18") {
    fl.open("../input_files/mcp2025b5e18_validation.txt");
    std::cout<<"filse in use are from ../input_files/mcp2025b5e18_validation.txt"<<std::endl;
  } else {
    std::cerr << "Failed to open the file for " << cintyp << "\n";
    exit(1);
  }
  return fl;
}

TString getOutputNamedQdx_official(const std::string& cintyp) {
  if (cintyp == "doData1e20") {
    return "../output_files/dQdx_data1e20_official.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/dQdx_mcp2025b5e18_official.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

TString getOutputNamedQdx_validation(const std::string& cintyp) {
  if (cintyp == "doData1e20") {
    return "../output_files/dQdx_data1e20_validation.root";
  } else if (cintyp == "doMC2025B5e18") {
    return "../output_files/dQdx_mcp2025b5e18_validation.root";
  } else {
    std::cerr << "Invalid cintyp provided. Exiting.\n";
    exit(1);
  }
}

