model ASF_Network_Simulation

import "./variables.gaml"
import "./map.gaml"
import "./farm.gaml"
import "./edge.gaml"
import "./network_functions.gaml"
import "./disease_transmission.gaml"
import "./culling_functions.gaml"
import "./statistics_functions.gaml"
import "./output_functions.gaml"
import "./gui_experiment.gaml"
import "./batch_experiment.gaml"

global {
    init {
        do load_spatial_data;
        do create_farms_batch;
        do set_connectivity_scenario;
        do initialize_network;
        do select_index_case;
        do output_initial_state; 
    }
    
    reflex weekly_update when: current_week < max_simulation_weeks {
        current_week <- current_week + 1;

        do update_statistics(); 
            
        do update_network_stergm();

        list<farm> newly_infected <- [];
        if num_infected_farms > 0 {
            newly_infected <- do_transmission();
            if !empty(newly_infected) {
                do process_infections(newly_infected);
            }
        }

        bool should_cull <- culling_enabled and 
                           first_infection_week >= 0 and 
                           (current_week - first_infection_week) >= culling_timing and
                           num_infected_farms > 0;
                           
        if should_cull {
            do apply_culling();
        }

        do calculate_stats();
        
        if num_infected_farms > peak_infected {
            peak_infected <- num_infected_farms;
            week_of_peak <- current_week;
        }

        int current_edges <- length(edge_connection);
        float current_mean_degree <- !empty(farm where (each.status != "removed")) ? 
              mean((farm where (each.status != "removed")) collect float(length(each.trading_partners))) : 0.0;
              
        if current_week <= 1 {
            write "";
            write "========================================";
            write "=== RUNNING SIMULATION ===";
            write "========================================";
        }
            
        write "Week " + current_week + " | Susceptible:" + num_susceptible_farms + 
              " Infected:" + num_infected_farms + " Remove:" + num_removed_farms +
              " | New Infected:" + length(newly_infected);
              
        write "Isolated farms: " + isolated + " (" + round(isolated * 100.0 / farms_count) + "%)";
              
//        write total_edges;
        
        if current_week >= max_simulation_weeks {
            do finalize();
        }
        
        if culling_enabled {
            if current_week >= culling_timing {
                do finalize();
            }
        }
    }
    
    reflex save_results when: (culling_enabled ? current_week >= culling_timing : current_week = max_simulation_weeks) {
        int total_infected_overall;
        int total_infected_small;
        int total_infected_medium;
        int total_infected_large;
        
    	if culling_enabled {
        	total_infected_overall <- num_removed_farms;
        	total_infected_small <- culled_by_type[1];
        	total_infected_medium <- culled_by_type[2];
        	total_infected_large <- culled_by_type[3];
        }
        else {
        	total_infected_overall <- length(farm where (each.infection_time >= 0));
        	total_infected_small <- length(farm where (each.infection_time >= 0 and each.farm_type = 1));
        	total_infected_medium <- length(farm where (each.infection_time >= 0 and each.farm_type = 2));
        	total_infected_large <- length(farm where (each.infection_time >= 0 and each.farm_type = 3));
        }
        
        
        bool csv_initialized <- file_exists(marker_file_path);
        
        if (!csv_initialized) {
            write "Initializing CSV file with headers";
            string header <- "cycle,connectivity_scenario,overall_mean_degree,total_infected_overall," +
                           "total_infected_small,total_infected_medium,total_infected_large," +
                           "peak_infected,week_of_peak";
            save header to: csv_file_path format: "csv" rewrite: true;
            
            save "initialized" to: marker_file_path format: "csv" rewrite: true;
        }
        
        string data_row <- "" + cycle + "," + 
                          connectivity_scenario + "," + 
                          overall_mean_degree + "," + 
                          total_infected_overall + "," + 
                          total_infected_small + "," + 
                          total_infected_medium + "," + 
                          total_infected_large + "," + 
                          peak_infected + "," + 
                          week_of_peak;
        
        save data_row to: csv_file_path format: "csv" rewrite: false;
    }
}