model ASF_Network_Simulation


global {
    file boundary_file <- file("../includes/data/shapefiles/hai_duong_boundary.shp");
    file roads_file <- file("../includes/data/shapefiles/hai_duong_roads.shp");
    file water_file <- file("../includes/data/shapefiles/hai_duong_water.shp");
    file waterways_file <- file("../includes/data/shapefiles/hai_duong_waterways.shp");
//    file farm_file <- file("../includes/data/farm_data/generated_farms_test.csv");
	file farm_file <- file("../includes/data/farm_data/farm.csv");
    
    geometry shape <- envelope(boundary_file + roads_file + water_file + waterways_file);
    
    float step <- 1 #week;
    int max_simulation_weeks <- 52;
    int current_week <- 0;
    int first_infection_week <- -1;
    
    float overall_mean_degree;
    int target_edges;
    
    bool edge_duration_enabled <- false;
    int edge_duration <- 26; // weeks
    float departure_rate <- 0.5;
    float homophily_coefficient <- 0.5;
    
    float transmission_probability <- 0.6;
    
    map<string, float> contact_rates <- [
        "1-1"::0.241, "1-2"::0.169, "1-3"::0.0,
        "2-1"::0.169, "2-2"::0.236, "2-3"::0.021,
        "3-1"::0.0,   "3-2"::0.021, "3-3"::0.021
    ];
    
    map<int, int> farms_by_type <- [1::0, 2::0, 3::0];
    map<int, int> infected_by_type <- [1::0, 2::0, 3::0];
    map<int, int> culled_by_type <- [1::0, 2::0, 3::0];
    
    string connectivity_scenario <- "baseline";
    bool culling_enabled <- false;
    int culling_timing <- 16; // weeks
    
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
    list<int> weekly_infected <- [];
    
    float mean_degree;
    int total_edges;
    int peak_infected <- 0;
    int week_of_peak <- 0;

    string csv_file_path <- "/Users/scott/DATN/ASF_ISP/asf_isp/results/scenario_connectivity_results.csv";
    string marker_file_path <- "/Users/scott/DATN/ASF_ISP/asf_isp/results/.csv_initialized";
    
    init {
        do load_spatial_data;
        do create_farms_batch;
        do set_connectivity_scenario;
        do initialize_network;
        do select_index_case;
        do output_initial_state; 
    }
    
    action load_spatial_data {
        create boundary from: boundary_file;
        create road from: roads_file;
        create water from: water_file;
        create waterways from: waterways_file;
    }
    
    action create_farms_batch {        
        create farm from: farm_file with: [
            location::point(float(get("x_coord")), float(get("y_coord"))),
            farm_id::string(get("farm_id")),
            farm_class::string(get("farm_class"))
        ] {
            if farm_class = "small" {
                farm_type <- 1;
            } else if farm_class = "medium" {
                farm_type <- 2;
            } else {
                farm_type <- 3;
            }
            
            status <- "susceptible";
            trading_partners <- [];
            infection_time <- -1;
            
            farms_by_type[farm_type] <- farms_by_type[farm_type] + 1;
        }
        
        num_susceptible_farms <- length(farm);
    }
    
    action set_connectivity_scenario {
        float scenario_value;
        if connectivity_scenario = "low" {
            scenario_value <- 0.5;
        } else if connectivity_scenario = "high" {
            scenario_value <- 1.0;
        } else {
            scenario_value <- 0.75; // baseline
        }
        
        overall_mean_degree <- scenario_value;
        target_edges <- int((length(farm) * overall_mean_degree) / 2.0);
    }
    
    action initialize_network {
        list<farm> all_farms <- list(farm);
        farms_count <- length(all_farms);
        
        int edges_created <- 0;
        int batch_size <- 5000;
        int max_total_attempts <- target_edges * 100;
        int total_attempts <- 0;
        
        list<pair<farm, farm>> potential_pairs <- [];
        
        int pair_generation_limit <- min(50000, farms_count * 20);
        int pairs_generated <- 0;
        
        loop while: pairs_generated < pair_generation_limit {
            farm f1 <- one_of(all_farms);
            farm f2 <- one_of(all_farms);
            
            if f1 != f2 and !(f1.trading_partners contains f2) {
                potential_pairs << (f1::f2);
                pairs_generated <- pairs_generated + 1;
            }
        }
        
        loop farm_pair over: potential_pairs {
            if edges_created >= target_edges { break; }
            
            farm f1 <- farm_pair.key;
            farm f2 <- farm_pair.value;
            
            if f1.trading_partners contains f2 { continue; }
            
            if length(f1.trading_partners) >= 100 or length(f2.trading_partners) >= 100 { 
                continue; 
            }
            
            string key <- string(f1.farm_type) + "-" + string(f2.farm_type);
            float base_prob <- contact_rates contains_key key ? contact_rates[key] : 0.0;
            
            float final_prob;
            if base_prob = 0.0 {
                final_prob <- 0.3;
            } else {
                final_prob <- min(0.98, base_prob * 10.0);
            }
            
            if f1.farm_type = f2.farm_type {
                final_prob <- min(0.99, final_prob * 1.5);
            }
            
            if flip(final_prob) {
                create edge_connection {
                    farm1 <- f1;
                    farm2 <- f2;
                    formation_time <- 0;
                    expected_duration <- edge_duration;
                }
                
                f1.trading_partners << f2;
                f2.trading_partners << f1;
                edges_created <- edges_created + 1;
            }
            
            total_attempts <- total_attempts + 1;
        }
        
		total_edges <- edges_created;
        
        if !empty(all_farms) {
            float actual_mean_degree <- mean(all_farms collect float(length(each.trading_partners)));
            mean_degree <- actual_mean_degree;
        }
        
        isolated <- length(all_farms where (length(each.trading_partners) = 0));
    }
    
    map<string, int> compute_type_pair_targets(map<int, list<farm>> farms_by_type_indexed) {
        map<string, int> targets;
        
        int n1 <- length(farms_by_type_indexed[1]);  // small
        int n2 <- length(farms_by_type_indexed[2]);  // medium  
        int n3 <- length(farms_by_type_indexed[3]);  // large
        
        int potential_11 <- int(n1 * (n1 - 1) / 2);  // within small
        int potential_22 <- int(n2 * (n2 - 1) / 2);  // within medium
        int potential_33 <- int(n3 * (n3 - 1) / 2);  // within large
        int potential_12 <- n1 * n2;                 // small-medium
        int potential_13 <- n1 * n3;                 // small-large
        int potential_23 <- n2 * n3;                 // medium-large
        
        float weight_11 <- contact_rates["1-1"] * (1.0 + homophily_coefficient);
        float weight_22 <- contact_rates["2-2"] * (1.0 + homophily_coefficient);
        float weight_33 <- contact_rates["3-3"] * (1.0 + homophily_coefficient);
        float weight_12 <- contact_rates["1-2"];
        float weight_13 <- contact_rates["1-3"];
        float weight_23 <- contact_rates["2-3"];
        
        float expected_11 <- potential_11 * weight_11;
        float expected_22 <- potential_22 * weight_22;
        float expected_33 <- potential_33 * weight_33;
        float expected_12 <- potential_12 * weight_12;
        float expected_13 <- potential_13 * weight_13;
        float expected_23 <- potential_23 * weight_23;
        
        float total_expected <- expected_11 + expected_22 + expected_33 + expected_12 + expected_13 + expected_23;
        
        if total_expected > 0 {
            float scale_factor <- target_edges / total_expected;
            
            targets["1-1"] <- int(expected_11 * scale_factor);
            targets["2-2"] <- int(expected_22 * scale_factor);
            targets["3-3"] <- int(expected_33 * scale_factor);
            targets["1-2"] <- int(expected_12 * scale_factor);
            targets["1-3"] <- int(expected_13 * scale_factor);
            targets["2-3"] <- int(expected_23 * scale_factor);
        } else {
            targets["1-1"] <- target_edges / 6;
            targets["2-2"] <- target_edges / 6;
            targets["3-3"] <- target_edges / 6;
            targets["1-2"] <- target_edges / 6;
            targets["1-3"] <- target_edges / 6;
            targets["2-3"] <- target_edges / 6;
        }
        
        return targets;
    }
    
    int create_edges_for_type_pair(string type_pair, int target_count, map<int, list<farm>> farms_by_type_indexed) {
        int type1;
        int type2;
        
        if type_pair = "1-1" {
            type1 <- 1; type2 <- 1;
        } else if type_pair = "1-2" {
            type1 <- 1; type2 <- 2;
        } else if type_pair = "1-3" {
            type1 <- 1; type2 <- 3;
        } else if type_pair = "2-2" {
            type1 <- 2; type2 <- 2;
        } else if type_pair = "2-3" {
            type1 <- 2; type2 <- 3;
        } else if type_pair = "3-3" {
            type1 <- 3; type2 <- 3;
        } else {
            type1 <- 1; type2 <- 1;
        }
        
        list<farm> farms1 <- copy(farms_by_type_indexed[type1]);
        list<farm> farms2;
        
        bool same_type <- (type1 = type2);
        if same_type {
            farms2 <- copy(farms1);
        } else {
            farms2 <- copy(farms_by_type_indexed[type2]);
        }
        
        if empty(farms1) or empty(farms2) {
            return 0;
        }
        
        farms1 <- shuffle(farms1);
        farms2 <- shuffle(farms2);
        
        int created <- 0;
        int max_attempts <- min(target_count * 5, length(farms1) * length(farms2));
        int attempts <- 0;
        
        int i1 <- 0;
        int i2 <- 0;
        
        loop while: created < target_count and attempts < max_attempts {
            attempts <- attempts + 1;
            
            farm f1 <- farms1[i1];
            farm f2 <- farms2[i2];
            
            if f1 = f2 {
                i2 <- (i2 + 1) mod length(farms2);
                if i2 = 0 { i1 <- (i1 + 1) mod length(farms1); }
                continue;
            }
            
            if f1.trading_partners contains f2 {
                i2 <- (i2 + 1) mod length(farms2);
                if i2 = 0 { i1 <- (i1 + 1) mod length(farms1); }
                continue;
            }
            
            if length(f1.trading_partners) >= 20 or length(f2.trading_partners) >= 20 {
                i2 <- (i2 + 1) mod length(farms2);
                if i2 = 0 { i1 <- (i1 + 1) mod length(farms1); }
                continue;
            }
            
            if flip(0.9) {
                create edge_connection {
                    farm1 <- f1;
                    farm2 <- f2;
                    formation_time <- 0;
                    expected_duration <- edge_duration;
                }
                
                f1.trading_partners << f2;
                f2.trading_partners << f1;
                created <- created + 1;
            }
            
            i2 <- (i2 + 1) mod length(farms2);
            if i2 = 0 { 
                i1 <- (i1 + 1) mod length(farms1); 
                if i1 mod 100 = 0 and !empty(farms2) {
                    farms2 <- shuffle(farms2);
                }
            }
        }
        
        return created;
    }
    
    action select_index_case {
        list<farm> medium_farms <- farm where (each.farm_type = 2);
        
        if !empty(medium_farms) {
            farm index_farm <- one_of(medium_farms);
            index_farm.status <- "infected";
            index_farm.infection_time <- 0;
            
            num_infected_farms <- 1;
            num_susceptible_farms <- num_susceptible_farms - 1;
            infected_by_type[2] <- 1;
            first_infection_week <- 0;
            
            index_id <- index_farm.farm_id;
        } else {
            error "No medium farms available for index case";
        }
    }

    reflex weekly_update when: current_week < max_simulation_weeks {
        current_week <- current_week + 1;

        list<farm> newly_infected <- [];
        if num_infected_farms > 0 {
            newly_infected <- do_transmission();
            if !empty(newly_infected) {
                do process_infections(newly_infected);
            }
        }

        do update_network_stergm();

        bool should_cull <- culling_enabled and 
                           first_infection_week >= 0 and 
                           (current_week - first_infection_week) >= culling_timing and
                           num_infected_farms > 0;
                           
        if should_cull {
            do apply_culling();
        }

        do calculate_stats();

        weekly_infected << num_infected_farms;
        
        if num_infected_farms > peak_infected {
            peak_infected <- num_infected_farms;
            week_of_peak <- current_week;
        }

        int current_edges <- length(edge_connection);
        float current_mean_degree <- !empty(farm where (each.status != "removed")) ? 
              mean((farm where (each.status != "removed")) collect float(length(each.trading_partners))) : 0.0;
            
        write "Week " + current_week + " | Susceptible:" + num_susceptible_farms + 
              " Infected:" + num_infected_farms + " Remove:" + num_removed_farms +
              " | New Infected:" + length(newly_infected);
        
        if current_week >= max_simulation_weeks {
            do finalize();
        }
    }
    
    list<farm> do_transmission {
    	list<farm> newly_infected <- [];
    	list<farm> infected_list <- farm where (each.status = "infected");
    
    	ask infected_list {
        	loop partner over: trading_partners {
            	if partner.status = "susceptible" and flip(transmission_probability) {
                	newly_infected << partner;
            	}
        	}
    	}
    
    	return remove_duplicates(newly_infected);
	}
    
    action process_infections(list<farm> newly_infected) {
        int count <- length(newly_infected);
        
        ask newly_infected {
            status <- "infected";
            infection_time <- current_week;
            infected_by_type[farm_type] <- infected_by_type[farm_type] + 1;
        }
        
        num_infected_farms <- num_infected_farms + count;
        num_susceptible_farms <- num_susceptible_farms - count;
        
        if first_infection_week < 0 {
            first_infection_week <- current_week;
        }
    }
    
    action update_network_stergm {
        list<edge_connection> edges_to_remove;
        
        if edge_duration_enabled {
            edges_to_remove <- edge_connection where 
                ((current_week - each.formation_time) >= each.expected_duration);
        } else {
            edges_to_remove <- edge_connection where flip(departure_rate);
        }
        
        int dissolved_count <- length(edges_to_remove);
        
        if dissolved_count > 0 {
            ask edges_to_remove {
                if farm1 != nil and farm2 != nil {
                    farm1.trading_partners >- farm2;
                    farm2.trading_partners >- farm1;
                }
                do die;
            }
        }
        
        int current_edges <- length(edge_connection);
        int edges_needed <- target_edges - current_edges;
        
        if edges_needed > 0 {
            int formed <- form_edges_aggressive(edges_needed);
        }
    }
    
    int form_edges_aggressive(int target_new_edges) {
        list<farm> active_farms <- farm where (each.status != "removed");
        if length(active_farms) < 2 { return 0; }
        
        int formed <- 0;
        int max_attempts <- target_new_edges * 20;
        int attempts <- 0;
        
        list<farm> available_farms <- active_farms where (length(each.trading_partners) < 50);
        
        if length(available_farms) < 2 {
            write "  WARNING: Most farms at degree limit, using all active farms";
            available_farms <- active_farms;
        }
        
        loop while: formed < target_new_edges and attempts < max_attempts {
            attempts <- attempts + 1;
            
            farm f1 <- one_of(available_farms);
            farm f2 <- one_of(available_farms);
            
            if f1 = f2 { continue; }
            if f1.trading_partners contains f2 { continue; }
            
            if length(f1.trading_partners) >= 100 or length(f2.trading_partners) >= 100 { 
                continue; 
            }
            
            string key <- string(f1.farm_type) + "-" + string(f2.farm_type);
            float base_prob <- contact_rates contains_key key ? contact_rates[key] : 0.0;
            
            float stergm_prob;
            if base_prob = 0.0 {
                stergm_prob <- 0.4;
            } else {
                stergm_prob <- min(0.95, base_prob * 8.0);
            }
            
            if f1.farm_type = f2.farm_type {
                stergm_prob <- min(0.98, stergm_prob * (1.0 + homophily_coefficient));
            }
            
            if flip(stergm_prob) {
                create edge_connection {
                    farm1 <- f1;
                    farm2 <- f2;
                    formation_time <- current_week;
                    expected_duration <- edge_duration;
                }
                
                f1.trading_partners << f2;
                f2.trading_partners << f1;
                formed <- formed + 1;
            }
        }
        return formed;
    }
    
    action apply_culling {
        list<farm> infected_list <- farm where (each.status = "infected");
        int num_culled <- length(infected_list);
        
        if num_culled = 0 { return; }
        
        write "=== CULLING INTERVENTION ===";
        write "Week " + current_week + ": Culling " + num_culled + " infected farms";
        
        ask infected_list {
            culled_by_type[farm_type] <- culled_by_type[farm_type] + 1;
            infected_by_type[farm_type] <- max(0, infected_by_type[farm_type] - 1);
        }

        list<edge_connection> edges_to_remove <- edge_connection where 
            (infected_list contains each.farm1 or infected_list contains each.farm2);
        
        ask edges_to_remove { do die; }

        ask infected_list {
            status <- "removed";
            ask trading_partners { trading_partners >- myself; }
            trading_partners <- [];
        }

        num_removed_farms <- num_removed_farms + num_culled;
        num_infected_farms <- 0;
        
        write "Culling complete. By type - Small: " + culled_by_type[1] + 
              ", Medium: " + culled_by_type[2] + 
              ", Large: " + culled_by_type[3];
    }
    
    
    action calculate_stats {
        
        list<farm> active <- farm where (each.status != "removed");
        
        if !empty(active) {
            list<int> degrees <- active collect length(each.trading_partners);
            mean_degree <- mean(degrees);
        } else {
            mean_degree <- 0.0;
        }
        
        total_edges <- length(edge_connection);
    }
    
    
    action output_initial_state {
        do calculate_stats();
        
        write "========================================";
        write "=== INITIAL STATE SUMMARY ===";
        write "========================================";
        write "Total farms: " + length(farm) + " (Small:" + farms_by_type[1] + 
              " Medium:" + farms_by_type[2] + " Large:" + farms_by_type[3] + ")";
        write "Network: " + total_edges + " edges, mean degree: " + 
              overall_mean_degree;
        write "Isolated farms: " + isolated + " (" + round(isolated * 100.0 / farms_count) + "%)";
        write "Index case selected: Medium farm " + index_id;
        write "Parameters: transmission=" + transmission_probability;
        write "Intervention: culling " + (culling_enabled ? "ON at week " + culling_timing : "OFF");
        
        if total_edges < (target_edges * 0.8) {
            write "WARNING: Network formation achieved only " + 
                  round(total_edges * 100.0 / target_edges) + "% of target connectivity";
        }
        
        write "";
        
        write "========================================";
        write "=== RUNNING SIMULATION ===";
        write "========================================";
    }
    
    action finalize {
        do calculate_stats();
        
        int total_affected <- num_removed_farms + num_infected_farms;
        int attack_rate <- round(total_affected * 100.0 / length(farm));
        
        write "";
        write "========================================";
        write "=== SIMULATION COMPLETE ===";
        write "========================================";
        write "Duration: " + current_week + " weeks";
        write "";
        write "FINAL STATUS:";
        write "  Susceptible: " + num_susceptible_farms + " (" + 
              round(num_susceptible_farms * 100.0 / length(farm)) + "%)";
        write "  Infected: " + num_infected_farms + " (" + 
              round(num_infected_farms * 100.0 / length(farm)) + "%)";
        write "  Removed: " + num_removed_farms + " (" + 
              round(num_removed_farms * 100.0 / length(farm)) + "%)";
        write "";
        write "EPIDEMIC SUMMARY:";
        write "  Peak infection: " + peak_infected + " farms at week " + week_of_peak;
        write "  Attack rate: " + attack_rate + "%";
        write "";
        write "INFECTIONS BY TYPE:";
        write "  Small: " + infected_by_type[1] + "/" + farms_by_type[1] + 
              " (" + round(infected_by_type[1] * 100.0 / farms_by_type[1]) + "%)";
        write "  Medium: " + infected_by_type[2] + "/" + farms_by_type[2] + 
              " (" + round(infected_by_type[2] * 100.0 / farms_by_type[2]) + "%)";
        write "  Large: " + infected_by_type[3] + "/" + farms_by_type[3] + 
              " (" + round(infected_by_type[3] * 100.0 / max(1, farms_by_type[3])) + "%)";
        
        if culling_enabled {
            int total_culled <- culled_by_type[1] + culled_by_type[2] + culled_by_type[3];
            write "";
            write "CULLING SUMMARY:";
            write "  Total culled: " + total_culled + " farms";
            write "  Small culled: " + culled_by_type[1];
            write "  Medium culled: " + culled_by_type[2];
            write "  Large culled: " + culled_by_type[3];
        }
        
        do pause;
    }
    
    action update_statistics {
        num_infected_farms <- length(farm where (each.status = "infected"));
        num_susceptible_farms <- length(farm where (each.status = "susceptible"));
        num_removed_farms <- length(farm where (each.status = "removed"));
        
        num_small_infected <- infected_by_type[1];
        num_medium_infected <- infected_by_type[2];
        num_large_infected <- infected_by_type[3];
        
        weekly_infected << num_infected_farms;
    }
    
    reflex save_results when: (current_week = max_simulation_weeks) {
        int total_infected_overall <- length(farm where (each.infection_time >= 0));
        int total_infected_small <- length(farm where (each.infection_time >= 0 and each.farm_type = 1));
        int total_infected_medium <- length(farm where (each.infection_time >= 0 and each.farm_type = 2));
        int total_infected_large <- length(farm where (each.infection_time >= 0 and each.farm_type = 3));
        
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

species boundary {
    aspect default {
        draw shape color: rgb(207, 207, 207) border: rgb(37, 37, 37) width: 1;
    }
}

species road {
    aspect default {
        draw shape color: rgb(60, 60, 60) width: 0.5;
    }
}

species water {
    aspect default {
        draw shape color: rgb(120, 172, 206) border: rgb(120, 172, 206);
    }
}

species waterways {
    aspect default {
        draw shape color: rgb(99, 140, 201) border: rgb(99, 140, 201);
    }
}

species farm {
    string farm_id;
    int farm_type; // 1=small, 2=medium, 3=large
    string farm_class;
    string status <- "susceptible";
    int infection_time <- -1;
    list<farm> trading_partners <- [];
    
    aspect default {
        float size;
        switch farm_type {
            match 1 { size <- 80.0; }
            match 2 { size <- 120.0; }
            match 3 { size <- 160.0; }
        }
        
        rgb color_farm;
        switch status {
            match "susceptible" { color_farm <- rgb(27, 158, 119); }
            match "infected" { color_farm <- rgb(255, 0, 255); }
            match "removed" { color_farm <- rgb(255, 0, 0); }
        }
        
        draw circle(size) color: color_farm border: #black width: 0.1;
    }
    
    aspect network_view {
        float size;
        switch farm_type {
            match 1 { size <- 80.0; }
            match 2 { size <- 120.0; }
            match 3 { size <- 160.0; }
        }
        
        rgb color_farm;
        switch status {
            match "susceptible" { color_farm <- rgb(27, 158, 119); }
            match "infected" { color_farm <- rgb(255, 0, 255); }
            match "removed" { color_farm <- rgb(255, 0, 0); }
        }
        
        draw circle(size) color: color_farm border: #black;
        
        loop partner over: trading_partners {
            draw line([location, partner.location]) color: rgb(100, 100, 100, 50) width: 1;
        }
    }
}

species edge_connection {
    farm farm1;
    farm farm2;
    int formation_time;
    int expected_duration;
}

experiment ASF_GUI_Simulation type: gui {
    parameter "Connectivity" var: connectivity_scenario among: ["low", "baseline", "high"];
    parameter "Enable Culling" var: culling_enabled;
    parameter "Culling Week" var: culling_timing min: 1 max: 52 step: 1;
    parameter "Enable Edge Duration" var: edge_duration_enabled;
    parameter "Edge Duration" var: edge_duration min: 6 max: 52 step: 1;
    parameter "Transmission Probability" var: transmission_probability min: 0.1 max: 1.0 step: 0.1;
    
    output synchronized: true {
        display map type: opengl 
        background: rgb(245, 247, 250)
        antialias: true {

            species boundary aspect: default refresh: false;
            species water aspect: default transparency: 0.5 refresh: false;
            species waterways aspect: default transparency: 0.5 refresh: false;
            species road aspect: default refresh: false;
            species farm aspect: default;
//            species farm aspect: network_view;
            
            // 2D Overlay
			overlay position: {1000, 1000} size: {14000, culling_enabled ? 28000 : 24000} 
            background: #white
            border: rgb(226, 232, 240)
            rounded: true
            transparency: 0.0 {
            
            // Panel dimensions
            float w <- 14000.0;
            float h <- 20000.0;

            // Responsive sizing (% of panel)
            float mx <- w * 0.05;
            float my <- h * 0.025;
            float cw <- w * 0.9;

            // Alignment anchors
            float cx <- w * 0.5;
            float lx <- mx;
            float rx <- w - mx;

            // Colors
            rgb green <- rgb(27, 158, 119);
            rgb purple <- rgb(255, 0, 255);
            rgb red <- rgb(255, 0, 0);
            rgb orange <- rgb(251, 146, 60);
            rgb text1 <- rgb(30, 41, 59);
            rgb text2 <- rgb(100, 116, 139);
            rgb card_green <- rgb(34, 197, 94, 20);
            rgb card_red <- rgb(239, 68, 68, 20);
            rgb badge_bg <- rgb(255, 255, 255, 200);
            rgb badge_border <- rgb(200, 200, 220);

            // Font sizes
            int fh <- int(w * 0.0015);  // header
            int fb <- int(w * 0.001);    // body
            int fs <- int(w * 0.001);    // small

            // Font definitions
            font header_bold <- font("Arial", fh, #bold);
            font body_bold <- font("Arial", fb, #bold);
            font body_reg <- font("Arial", fb);
            font small_bold <- font("Arial", fs, #bold);
            font small_reg <- font("Arial", fs);

            // Layout metrics
            float sh <- h * 0.1;
            float ls <- h * 0.06;
            float y <- my;

            // Header
            draw "ASF MODEL STATUS" at: {cx, y + sh/2}
                color: #black font: header_bold anchor: #center;
            y <- y + sh + ls * 0.5;

            // Week badge
            float bw <- w * 0.3;
            float bh <- h * 0.05;
            draw "WEEK " + current_week at: {cx, edge_duration_enabled ? y + bh/10 : y + bh/2}
                color: text2 font: small_bold anchor: #center;
            
            // Edge duration
            if (edge_duration_enabled) {
            draw "Edge duration: " + edge_duration + " (week)"at: {cx, y + bh * 1.3}
                color: text2 font: small_bold anchor: #center;  
            }
            y <- y + bh + ls;
            

            // Total Farms
            float fh2 <- sh * 0.8;
            draw "TOTAL FARMS" at: {lx, y + fh2/2}
                color: text1 font: body_bold anchor: #left_center;
            draw string(length(farm)) at: {rx - mx, y + fh2/2}
                color: text1 font: header_bold anchor: #right_center;
            y <- y + fh2 + ls * 0.5;

            // First divider
            draw line([{lx, y}, {rx, y}]) color: #black width: 1;
            y <- y + ls;

            // Infection Status
            draw "INFECTION STATUS" at: {lx, y}
                color: text2 font: small_bold anchor: #left_center;
            y <- y + ls;

            // Status cards
            float cardw <- (cw - mx) * 0.45;
            float cardh <- sh * 1.2;
            float lcx <- lx + cardw/2;
            float rcx <- rx - cardw/2;
            float mcx <- rcx - lcx;
            float card_y <- y + cardh/3;
            float text_y1 <- y + cardh * 0.075;
            float text_y2 <- y + cardh * 0.56;

            // Susceptible card
            draw rectangle(cardw, cardh) at: {lcx, card_y}
                color: card_green border: green width: 2;
            draw "Susceptible" at: {lcx, text_y1}
                color: green font: body_reg anchor: #center;
            draw string(num_susceptible_farms) at: {lcx, text_y2}
                color: green font: header_bold anchor: #center;

            // Infected card
            draw rectangle(cardw, cardh) at: {rcx, card_y}
                color: card_red border: purple width: 2;
            draw "Infected" at: {rcx, text_y1}
                color: purple font: body_reg anchor: #center;
            draw string(num_infected_farms) at: {rcx, text_y2}
                color: purple font: header_bold anchor: #center;
                
            // Removed card
            if (culling_enabled) {
            	draw rectangle(cardw, cardh * 1.8) at: {mcx, card_y + cardh * 1.7}
                	color: card_red border: red width: 2;
            	draw "Removed" at: {mcx, text_y1 * 1.31}
                	color: red font: body_reg anchor: #center;
               	draw "(week: " + culling_timing + ")" at: {mcx, text_y1 * 1.4}
                	color: red font: small_reg anchor: #center;
            	draw string(num_removed_farms) at: {mcx, text_y2 * 1.4}
                	color: red font: header_bold anchor: #center;
            }

            y <- culling_enabled ? y + cardh + ls + cardh * 1.7 : y + cardh + ls;

            // Second divider
            draw line([{lx, y}, {rx, y}]) color: #black width: 1;
            y <- y + ls;

            // Breakdown by Size
            draw "BREAKDOWN BY SIZE" at: {lx, y}
                color: text2 font: small_bold anchor: #left_center;
            y <- y + ls;

            // Bar chart setup
            int maxv <- culling_enabled ? (current_week < culling_timing ? num_small_infected + num_medium_infected + num_large_infected : num_small_culled + num_medium_culled + num_large_culled) : num_small_infected + num_medium_infected + num_large_infected;
            float barh <- h * 0.01;
            float barw <- cw * 0.9;
            float bar_offset <- 700.0;

            // Small farms bar
            float small_w <- (culling_enabled ? (current_week < culling_timing ? num_small_infected / maxv : num_small_culled / maxv) : num_small_infected / maxv) * barw;
            draw "Small" at: {lx, y}
                color: text1 font: small_reg anchor: #left_center;
            draw rectangle(small_w, barh) at: {lx + small_w/2, y + bar_offset}
                color: culling_enabled ? (current_week < culling_timing ? purple : red ) : purple;
            draw string(culling_enabled ? (current_week < culling_timing ? num_small_infected : num_small_culled) : num_small_infected) at: {rx, y}
                color: text1 font: small_bold anchor: #right_center;
            y <- y + ls * 1.5;

            // Medium farms bar
            float medium_w <- (culling_enabled ? (current_week < culling_timing ? num_medium_infected / maxv : num_medium_culled / maxv) : num_medium_infected / maxv) * barw;
            draw "Medium" at: {lx, y}
                color: text1 font: small_reg anchor: #left_center;
            draw rectangle(medium_w, barh) at: {lx + medium_w/2, y + bar_offset}
                color: culling_enabled ? (current_week < culling_timing ? purple : red ) : purple;
            draw string(culling_enabled ? (current_week < culling_timing ? num_medium_infected : num_medium_culled) : num_medium_infected) at: {rx, y}
                color: text1 font: small_bold anchor: #right_center;
            y <- y + ls * 1.5;

            // Large farms bar
            float large_w <- (culling_enabled ? (current_week < culling_timing ? num_large_infected / maxv : num_large_culled / maxv) : num_large_infected / maxv) * barw;
            draw "Large" at: {lx, y}
                color: text1 font: small_reg anchor: #left_center;
            draw rectangle(large_w, barh) at: {lx + large_w/2, y + bar_offset}
                color: culling_enabled ? (current_week < culling_timing ? purple : red ) : purple;
            draw string(culling_enabled ? (current_week < culling_timing ? num_large_infected : num_large_culled) : num_large_infected) at: {rx, y}
                color: text1 font: small_bold anchor: #right_center;
            y <- y + ls * 1.2;

            // Third divider
            draw line([{lx, y}, {rx, y}]) color: #black width: 1;
            y <- y + ls * 0.5;

            // Scenario section
            float sch <- sh * 0.7;
            draw "SCENARIO" at: {lx, y + sch * 0.3}
                color: text2 font: small_bold anchor: #left_center;

            string scenario_text <- connectivity_scenario + " (" + string(round(mean_degree * 100) / 100) + ")";
            draw scenario_text at: {lx, y + sch}
                color: text1 font: body_bold anchor: #left_center;

            // Infection rate badge
            if (length(farm) > 0) {
            	int total_culled_farms <- num_small_culled + num_medium_culled +num_large_culled;
                float rate <- culling_enabled ? (total_culled_farms / length(farm)) * 100 : (num_infected_farms / length(farm)) * 100;
                rgb rate_color <- orange;
                float radius <- w * 0.08;
                float badge_x <- rx - radius - mx;
                float badge_y <- y + sch/1.5;

                draw circle(radius) at: {badge_x, badge_y}
                    color: #white border: rate_color width: 2;
                draw string(round(rate)) + "%" at: {badge_x, badge_y}
                    color: rate_color font: small_bold anchor: #center;
            }
            }
        }
        
        monitor "Week" value: current_week;
        monitor "Susceptible" value: num_susceptible_farms;
        monitor "Infected" value: num_infected_farms;
        monitor "Removed" value: num_removed_farms;
        monitor "Mean Degree" value: round(mean_degree * 100) / 100;
        monitor "Current Edges" value: total_edges;
        monitor "Small Infected" value: num_small_infected;
        monitor "Medium Infected" value: num_medium_infected;
        monitor "Large Infected" value: num_large_infected;
    }
}


experiment Batch_Connectivity_Scenarios type: batch repeat: 200 keep_seed: false parallel: 5 until: (current_week >= max_simulation_weeks) {
    
    parameter "Enable Culling" var: culling_enabled init: false;
    parameter "Culling Week" var: culling_timing init: 16;
    parameter "Enable Edge Duration" var: edge_duration_enabled init: false;
    parameter "Edge Duration" var: edge_duration init: 26;
    parameter "Transmission Probability" var: transmission_probability init: 0.6;
    parameter "Connectivity" var: connectivity_scenario init: "baseline";
}