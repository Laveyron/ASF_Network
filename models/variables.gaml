model ASF_Network_Simulation

global {
    file boundary_file <- file("../includes/data/shapefiles/hai_duong_boundary.shp");
    file roads_file <- file("../includes/data/shapefiles/hai_duong_roads.shp");
    file water_file <- file("../includes/data/shapefiles/hai_duong_water.shp");
    file waterways_file <- file("../includes/data/shapefiles/hai_duong_waterways.shp");
    file farm_file <- file("../includes/data/farm_data/farm.csv");
//    file farm_file <- file("../includes/data/farm_data/generated_farms_test.csv"); 
    
    geometry shape <- envelope(boundary_file + roads_file + water_file + waterways_file);
    
    float step <- 1 #week;
    int max_simulation_weeks <- 52;
    int current_week <- 0;
    int first_infection_week <- -1;
    
    float overall_mean_degree;
    int target_edges;
    
    bool edge_duration_enabled <- false;
    int edge_duration <- 26;
    float dissolution_rate <- 1.0;

	map<string, float> contact_rates <- [
        "1-1"::0.241, "1-2"::0.169, "1-3"::0.0,
        "2-1"::0.241, "2-2"::0.236, "2-3"::0.021,
        "3-1"::0.0,   "3-2"::0.236, "3-3"::0.021
    ];
    
    map<int, int> farms_by_type <- [1::0, 2::0, 3::0];
    map<int, int> infected_by_type <- [1::0, 2::0, 3::0];
    map<int, int> culled_by_type <- [1::0, 2::0, 3::0];
    
    string connectivity_scenario <- "Cơ sở";
    bool culling_enabled <- false;
    int culling_timing <- 16;
    
    float k0 <- 1.0;           
    float r0 <- 3.0;      
    float alpha <- 2.27;

    float transmission_probability <- 0.6;
    
    int farms_count;
    int isolated;
    string index_id;
    int num_infected_farms;
    int num_susceptible_farms;
    int num_removed_farms;
    int num_small_infected;
    int num_medium_infected;
    int num_large_infected;
    int num_small_culled;
    int num_medium_culled;
    int num_large_culled;
    
    float mean_degree;
    int total_edges;
    int peak_infected <- 0;
    int week_of_peak <- 0;
	
	string experiment_id;
    string csv_file_path <- "../results/kernel_baseline_results.csv";
    string marker_file_path <- "../results/.csv_initialized";
}
