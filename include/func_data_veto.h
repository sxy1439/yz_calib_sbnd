
// plane 0

bool InVeto_region_eastTPC_I0(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 230, 280},      // Middle-bottom horizontal region (Region 2)
    {-195, -190, 250, 255},
    {-200, -195, 345, 455},      // Small bottom-right extension (Region 3)
    {-200, -185, 495, 500},       // Thin vertical strip near the far-right side (Region 4)
    {-10, 10, 495, 500},        // Large vertical block near the top corner (Region 5)
    {190, 200, 460, 500},
    {195, 200, 235, 285},        // Horizontal region in the middle (Region 6)
    {190, 195, 255, 275},           // Vertical region on left side (Region 7)
    {185, 200, 0, 15},            // Region near bottom left (Region 8)       
    {-20, 25, 5, 10},
    {-20, 0, 10, 15},
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_I0(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -180, 0, 10},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 235, 270},      // Middle-bottom horizontal region (Region 2)
    {-195, -185, 250, 260},
    {-200, -195, 345, 490},       // Thin vertical strip near the far-right side (Region 3)
    {-5, 5, 495, 500},
    {185, 200, 480, 500},        // Large vertical block near the top corner (Region 4)
    {190, 200, 245, 260},
    {185, 200, 0, 15},           // Vertical region on left side (Region 6)
    {-25, 10, 5, 10},            // Region near bottom left (Region 8)
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


// plane 1

bool InVeto_region_eastTPC_I1(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 250, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -195, 425, 460},      // Small bottom-right extension (Region 3)
    {-200, -180, 480, 500},       // Thin vertical strip near the far-right side (Region 4)
    {180, 200, 460, 490},        // Large vertical block near the top corner (Region 5)
    {190, 200, 490, 500},
    {185, 200, 245, 270},        // Horizontal region in the middle (Region 6)
    {185, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-15, 10, 5, 10},            // Region near bottom left (Region 8)       
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_I1(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -180, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -185, 245, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -175, 495, 500},       // Thin vertical strip near the far-right side (Region 3)
    {185, 200, 485, 500},        // Large vertical block near the top corner (Region 4)
    //{185, 200, 235, 265},        // Horizontal region in the middle (Region 5)
    {190, 200, 250, 265},
    {185, 200, 0, 15},           // Vertical region on left side (Region 6)
    //{170, 190, 10, 15},
    //{-5, 5, 0, 15},            // Region near bottom left (Region 7)
    {-25, 15, 5, 10},            // Region near bottom left (Region 8)
    {-15, 5, 10, 15},
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}




// plane 2

bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 240, 260},      // Middle-bottom horizontal region (Region 2)
    //{-200, -195, 445, 455},      // Small bottom-right extension (Region 3)
    {-200, -175, 490, 500},       // Thin vertical strip near the far-right side (Region 4)
    {180, 200, 460, 490},        // Large vertical block near the top corner (Region 5)
    {190, 200, 490, 500},
    {185, 200, 230, 265},        // Horizontal region in the middle (Region 6)
    {180, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-30, 20, 5, 15},            // Region near bottom left (Region 8)       
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 235, 265},      // Middle-bottom horizontal region (Region 2)
    {-200, -175, 495, 500},       // Thin vertical strip near the far-right side (Region 3)
    {180, 200, 490, 500},        // Large vertical block near the top corner (Region 4)
    {185, 200, 235, 265},        // Horizontal region in the middle (Region 5)
    {190, 200, 235, 270},
    {190, 200, 0, 15},           // Vertical region on left side (Region 6)
    {170, 190, 10, 15},
    //{-5, 5, 0, 15},            // Region near bottom left (Region 7)
    {-30, 45, 5, 10},            // Region near bottom left (Region 8)
    {-10, 10, 10, 15},
    {-200, 200, 0, 5}
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}




bool InVeto_region_X(double x) {
  vector<vector<double>> vetoed_regions = {
    {-200, -195},      // anode 1
    {-5, 5},          // cathode
    {195, 200}       // anode 2
    
  };

  // see if the point falls in the veto region
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double x_min = region[0];
    double x_max = region[1];
    
    if (x >= x_min && x <= x_max) {
      return true;
    }
  } 
  return false; 
}

/*
// Test 1 (Subset ~20k tracks : presented on Jan 30, 2025)

// subset of presented one 
bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -190, 235, 265},      // Middle-bottom horizontal region (Region 2)
    {-200, -190, 440, 485},      // Small bottom-right extension (Region 3)
    {-200, -175, 490, 500},       // Thin vertical strip near the far-right side (Region 4)
    {180, 200, 400, 500},        // Large vertical block near the top corner (Region 5)
    {185, 200, 230, 270},        // Horizontal region in the middle (Region 6)
    {185, 200, 0, 20},           // Vertical region on left side (Region 7)
    {-60, 45, 0, 10},            // Region near bottom left (Region 8)
    {-10, 5, 10, 15}             // Region near bottom left (Region 8, set 2)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -190, 240, 280},      // Middle-bottom horizontal region (Region 2)
    {-200, -165, 485, 500},       // Thin vertical strip near the far-right side (Region 3)
    {180, 200, 475, 500},        // Large vertical block near the top corner (Region 4)
    {190, 200, 230, 275},        // Horizontal region in the middle (Region 5)
    {185, 200, 0, 20},           // Vertical region on left side (Region 6)
    {-60, 40, 0, 10},            // Region near bottom left (Region 8)
    {-10, 5, 10, 15}            // Region near bottom left (Region 8, set 2)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}
*/


/*
//presented on Jan 30, 2025

bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -180, 0, 20},         // Bottom-left horizontal region (Region 1)
    {-200, -190, 225, 270},      // Middle-bottom horizontal region (Region 2)
    {-200, -190, 420, 485},      // Small bottom-right extension (Region 3)
    {-200, -175, 490, 500},       // Thin vertical strip near the far-right side (Region 4)
    {170, 200, 310, 500},        // Large vertical block near the top corner (Region 5)
    {180, 200, 230, 270},        // Horizontal region in the middle (Region 6)
    {190, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-60, 40, 0, 10}            // Region near bottom left (Region 8)       
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -180, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -190, 240, 280},      // Middle-bottom horizontal region (Region 2)
    {-200, -130, 480, 500},       // Thin vertical strip near the far-right side (Region 3)
    {170, 200, 470, 500},        // Large vertical block near the top corner (Region 4)
    {190, 200, 240, 265},        // Horizontal region in the middle (Region 5)
    {190, 200, 0, 15},           // Vertical region on left side (Region 6)
    //{-5, 5, 0, 15}            // Region near bottom left (Region 7)
    {-60, 40, 0, 10}            // Region near bottom left (Region 8)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}
*/



// Test 2 (Subset ~21k tracks) : Channel to channel mapping applied


/*
bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 240, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -195, 445, 455},      // Small bottom-right extension (Region 3)
    {-200, -185, 495, 500},       // Thin vertical strip near the far-right side (Region 4)
    {190, 200, 460, 500},        // Large vertical block near the top corner (Region 5)
    {185, 200, 230, 265},        // Horizontal region in the middle (Region 6)
    {180, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-55, 45, 0, 10}            // Region near bottom left (Region 8)       
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 230, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -180, 485, 500},       // Thin vertical strip near the far-right side (Region 3)
    {190, 200, 490, 500},        // Large vertical block near the top corner (Region 4)
    {190, 200, 235, 275},        // Horizontal region in the middle (Region 5)
    {190, 200, 0, 15},           // Vertical region on left side (Region 6)
    //{-5, 5, 0, 15}            // Region near bottom left (Region 7)
    {-50, 35, 0, 10}            // Region near bottom left (Region 8)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}

*/


// Test 1 (full stats)

/*
bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 240, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -195, 420, 485},      // Small bottom-right extension (Region 3)
    {-200, -180, 495, 500},       // Thin vertical strip near the far-right side (Region 4)
    {185, 200, 460, 500},        // Large vertical block near the top corner (Region 5)
    {180, 200, 245, 265},        // Horizontal region in the middle (Region 6)
    {185, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-55, 45, 0, 10},            // Region near bottom left (Region 8)
    {-200, 200, 65, 70}             // Region of the dead channel (Region 9)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -190, 235, 270},      // Middle-bottom horizontal region (Region 2)
    {-200, -195, 420, 465},      // Small bottom-right extension (Region 3)
    {-200, -145, 480, 500},       // Thin vertical strip near the far-right side (Region 4)
    {175, 200, 475, 500},        // Large vertical block near the top corner (Region 5)
    {185, 200, 230, 270},        // Horizontal region in the middle (Region 6)
    {185, 200, 0, 15},           // Vertical region on left side (Region 7)
    {-30, 25, 0, 10}            // Region near bottom left (Region 8)       
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}
*/


// Test 3
/*
bool InVeto_region_eastTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 25, 500},      // Middle-bottom horizontal region (Region 2)
    //{-200, -195, 445, 455},      // Small bottom-right extension (Region 3)
    //{-200, -185, 495, 500},       // Thin vertical strip near the far-right side (Region 4)
    {185, 200, 450, 500},        // Large vertical block near the top corner (Region 5)
    {185, 200, 230, 260},        // Horizontal region in the middle (Region 6)
    {180, 200, 10, 20},           // Vertical region on left side (Region 7)
    {-185, 200, 0, 10}            // Region near bottom left (Region 8)       
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}


bool InVeto_region_westTPC_C(double y, double z) {
  // format considering {y_min, y_max, z_min, z_max}
  vector<vector<double>> vetoed_regions = {
    {-200, -185, 0, 15},         // Bottom-left horizontal region (Region 1)
    {-200, -195, 230, 260},      // Middle-bottom horizontal region (Region 2)
    {-200, -180, 485, 500},       // Thin vertical strip near the far-right side (Region 3)
    {190, 200, 490, 500},        // Large vertical block near the top corner (Region 4)
    {190, 200, 235, 275},        // Horizontal region in the middle (Region 5)
    {190, 200, 0, 15},           // Vertical region on left side (Region 6)
    //{-5, 5, 0, 15}            // Region near bottom left (Region 7)
    {-50, 35, 0, 10}            // Region near bottom left (Region 8)
  };

  // to see if the point falls in the veto region
  //for (const auto& region : vetoed_regions) { // using element in the container method
  for (size_t i = 0; i < vetoed_regions.size(); ++i) {
    vector<double> region = vetoed_regions[i];
    
    double y_min = region[0];
    double y_max = region[1];
    double z_min = region[2];
    double z_max = region[3];
    
    if (y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
      return true;
    }
  } 
  return false; 
}
*/
