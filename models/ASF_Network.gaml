model ASF_Network_Simulation

global {
    // File paths
    file boundary_file <- file("../includes/data/shapefiles/hai_duong_boundary.shp");
    file roads_file <- file("../includes/data/shapefiles/hai_duong_roads.shp");
    file water_file <- file("../includes/data/shapefiles/hai_duong_water.shp");
    file waterways_file <- file("../includes/data/shapefiles/hai_duong_waterways.shp");
//    file farm_file <- file("../includes/data/farm_data/farm.csv");
    file farm_file <- file("../includes/data/farm_data/generated_farms_test.csv");
    
    geometry shape <- envelope(boundary_file + roads_file + water_file + waterways_file);
    
    // Time parameters
    float step <- 1 #week;
    int max_simulation_weeks <- 52;
    int current_week <- 0;
    int first_infection_week <- -1;
    
    // STERGM Parameters
    float theta_formation_edges <- -3.0;      // Edge formation intercept
    float theta_formation_homophily <- 1.5;   // Same-type preference
    float theta_dissolution_edges <- -2.0;    // Edge persistence intercept  
    float theta_dissolution_duration <- 0.05; // Duration dependence
    
    // Edge duration parameters
    bool edge_duration_enabled <- false;      // Toggle for edge duration feature
    int edge_duration <- 26;              // Minimum weeks to sustain connection (default: 26)
    float departure_rate <- 0.5;              // Random dissolution rate when duration not enabled
    
    // Network statistics for monitoring
    float target_mean_degree <- 0.75;  // TARGET MEAN DEGREE - PRIMARY PARAMETER (baseline)
    float target_density;               // Calculated from mean degree
    int target_edges;                   // Calculated from mean degree
    float current_density;
    float current_mean_degree <- 0.0;
    
    // Fixed epidemiological parameters
    float transmission_probability <- 0.6;
    
    // Fixed contact rates matrix (normalized as probabilities)
    map<string, float> contact_rates <- [
        "1-1"::0.241, "1-2"::0.169, "1-3"::0.0,
        "2-1"::0.169, "2-2"::0.236, "2-3"::0.021,
        "3-1"::0.0,   "3-2"::0.021, "3-3"::0.021
    ];
    
    // Tracking variables
    map<int, int> farms_by_type <- [1::0, 2::0, 3::0];
    map<int, int> infected_by_type <- [1::0, 2::0, 3::0];
    map<int, int> culled_by_type <- [1::0, 2::0, 3::0];
    
    // Intervention parameters
    string connectivity_scenario <- "baseline";
    bool culling_enabled <- false;
    int culling_timing <- 16;
    
    // Statistics
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
    
    // STERGM-specific tracking
    map<string, float> network_statistics;
    list<float> density_history <- [];
    list<float> transitivity_history <- [];
    
    string csv_file_path <- "/Users/scott/DATN/ASF_ISP/asf_isp/results/stergm_results.csv";
    string marker_file_path <- "/Users/scott/DATN/ASF_ISP/asf_isp/results/.stergm_initialized";
    
    init {
        do load_spatial_data;
        do create_farms_batch;
        do set_stergm_parameters;
        do initialize_network_stergm;
        do select_index_case;
        do compute_network_statistics;
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
        farms_count <- length(farm);
    }
    
    action set_stergm_parameters {
        // Set target mean degree based on connectivity scenario - FIXED VALUES
        if connectivity_scenario = "low" {
            target_mean_degree <- 0.5;    // Fixed low connectivity
            theta_formation_edges <- -2.5;  // Less negative for better edge formation
            theta_dissolution_edges <- -3.0; // More stable edges
        } else if connectivity_scenario = "high" {
            target_mean_degree <- 1.0;     // Fixed high connectivity
            theta_formation_edges <- -2.0;  // Even less negative
            theta_dissolution_edges <- -3.5; // Very stable edges
        } else { // baseline
            target_mean_degree <- 0.75;    // Fixed baseline connectivity
            theta_formation_edges <- -2.3;  // Adjusted for better formation
            theta_dissolution_edges <- -3.2; // Stable edges
        }
        
        // Calculate target density and edges from mean degree
        int n <- length(farm);
        // density = mean_degree / (n - 1)
        target_density <- target_mean_degree / (n - 1.0);
        // total_edges = (n * mean_degree) / 2
        target_edges <- int((n * target_mean_degree) / 2.0);
        
        write "Target Mean Degree (FIXED): " + target_mean_degree;
        write "Calculated Target Density: " + target_density;
        write "Calculated Target Edges: " + target_edges;
    }
    
    action initialize_network_stergm {
        // Initialize network using STERGM formation model
        list<farm> all_farms <- list(farm);
        int n <- length(all_farms);
        
        write "Initializing STERGM network with target mean degree: " + target_mean_degree;
        
        // Create initial edges based on formation model
        int edges_created <- 0;
        int attempts <- 0;
        int max_attempts <- target_edges * 50; // More attempts for sparse networks
        
        // Keep trying until we reach at least 80% of target
        loop while: edges_created < target_edges * 0.8 and attempts < max_attempts {
            attempts <- attempts + 1;
            
            farm f1 <- one_of(all_farms);
            farm f2 <- one_of(all_farms);
            
            if f1 = f2 or f1.trading_partners contains f2 {
                continue;
            }
            
            // Calculate formation probability using STERGM model
            float log_odds <- calculate_formation_log_odds(f1, f2, true);
            float prob <- 1.0 / (1.0 + exp(-log_odds));
            
            // Boost probability during initialization to reach target
            prob <- min(0.95, prob * 2.0);
            
            if flip(prob) {
                create edge_connection {
                    farm1 <- f1;
                    farm2 <- f2;
                    formation_time <- 0;
                    edge_duration <- 0;
                }
                
                f1.trading_partners << f2;
                f2.trading_partners << f1;
                edges_created <- edges_created + 1;
            }
        }
        
        // If still below target, force additional random edges
        if edges_created < target_edges * 0.8 {
            write "  Forcing additional edges to reach target...";
            loop while: edges_created < target_edges {
                farm f1 <- one_of(all_farms);
                farm f2 <- one_of(all_farms);
                
                if f1 != f2 and !(f1.trading_partners contains f2) {
                    create edge_connection {
                        farm1 <- f1;
                        farm2 <- f2;
                        formation_time <- 0;
                        edge_duration <- 0;
                    }
                    
                    f1.trading_partners << f2;
                    f2.trading_partners << f1;
                    edges_created <- edges_created + 1;
                }
            }
        }
        
        total_edges <- edges_created;
        current_density <- (2.0 * total_edges) / (n * (n - 1));
        
        // Calculate actual mean degree
        list<int> degrees <- all_farms collect length(each.trading_partners);
        current_mean_degree <- mean(degrees);
        
        write "Initial network created: " + total_edges + " edges";
        write "Actual mean degree: " + current_mean_degree + " (target: " + target_mean_degree + ")";
        write "Achievement rate: " + round((current_mean_degree / target_mean_degree) * 100) + "%";
        
        // Calculate initial statistics
        do compute_network_statistics;
    }
    
    float calculate_formation_log_odds(farm f1, farm f2, bool is_initial <- false) {
        float log_odds <- theta_formation_edges;
        
        // Homophily effect
        if f1.farm_type = f2.farm_type {
            log_odds <- log_odds + theta_formation_homophily;
        }
        
        // Contact rate effect (from fixed matrix)
        string key <- string(min(f1.farm_type, f2.farm_type)) + "-" + 
                     string(max(f1.farm_type, f2.farm_type));
        float contact_weight <- contact_rates contains_key key ? contact_rates[key] : 0.1;
        
        // Stronger influence of contact rates for sparse networks
        log_odds <- log_odds + ln(max(0.01, contact_weight)) * 3.0;  // Increased from 2.0
        
        // Very light degree penalty for sparse networks
        if !is_initial {
            int current_degree_sum <- length(f1.trading_partners) + length(f2.trading_partners);
            // Even lighter penalty for sparse networks
            float degree_penalty <- -0.002 * current_degree_sum;  // Reduced from -0.005
            log_odds <- log_odds + degree_penalty;
        }
        
        return log_odds;
    }
    
    float calculate_dissolution_log_odds(edge_connection e) {
        // If edge duration is enabled and edge hasn't reached minimum duration, prevent dissolution
        if edge_duration_enabled {
            int current_duration <- current_week - e.formation_time;
            if current_duration < edge_duration {
                // Return very negative log odds to prevent dissolution
                return -10.0;  // Effectively probability ~ 0
            }
        }
        
        // Normal STERGM dissolution calculation
        float log_odds <- theta_dissolution_edges;
        
        // Duration dependence (edges become more stable over time)
        int duration <- current_week - e.formation_time;
        log_odds <- log_odds - theta_dissolution_duration * duration;
        
        // Homophily effect on dissolution (same types maintain connections longer)
        if e.farm1.farm_type = e.farm2.farm_type {
            log_odds <- log_odds - 0.5; // Less likely to dissolve
        }
        
        // Additional stability when below target mean degree
        if current_mean_degree < target_mean_degree * 0.9 {
            log_odds <- log_odds - 1.0; // Make edges more stable when below target
        }
        
        return log_odds;
    }
    
    reflex weekly_update when: current_week < max_simulation_weeks {
        current_week <- current_week + 1;
        
        // STERGM network evolution
        do evolve_network_stergm;
        
        // Disease transmission
        list<farm> newly_infected <- [];
        if num_infected_farms > 0 {
            newly_infected <- do_transmission;
            if !empty(newly_infected) {
                do process_infections(newly_infected);
            }
        }
        
        // Culling intervention
        bool should_cull <- culling_enabled and 
                           first_infection_week >= 0 and 
                           (current_week - first_infection_week) >= culling_timing and
                           num_infected_farms > 0;
                           
        if should_cull {
            do apply_culling;
        }
        
        // Update statistics
        do compute_network_statistics;
        weekly_infected << num_infected_farms;
        
        if num_infected_farms > peak_infected {
            peak_infected <- num_infected_farms;
            week_of_peak <- current_week;
        }
        
        // Log progress
        write "Week " + current_week + 
              " | S:" + num_susceptible_farms + 
              " I:" + num_infected_farms + 
              " R:" + num_removed_farms +
              " | Mean degree: " + round(current_mean_degree * 100) / 100 +
              " (target: " + target_mean_degree + ")" +
              " | Edges: " + total_edges + 
              " | New infections: " + length(newly_infected);
        
        if current_week >= max_simulation_weeks {
            do finalize;
        }
    }
    
    action evolve_network_stergm {
        // STERGM evolution: separate dissolution and formation processes
        
        // Phase 1: Edge Dissolution
        list<edge_connection> edges_to_evaluate <- list(edge_connection);
        list<edge_connection> edges_to_remove <- [];
        
        ask edges_to_evaluate {
            // Skip edges involving removed farms
            if farm1.status = "removed" or farm2.status = "removed" {
                edges_to_remove << self;
            } else {
                // Calculate dissolution probability
                float log_odds <- myself.calculate_dissolution_log_odds(self);
                float prob_dissolve <- 1.0 / (1.0 + exp(-log_odds));
                
                if flip(prob_dissolve) {
                    edges_to_remove << self;
                } else {
                    // Update edge duration
                    edge_duration <- edge_duration + 1;
                }
            }
        }
        
        // Remove dissolved edges
        ask edges_to_remove {
            if farm1 != nil and farm2 != nil {
                farm1.trading_partners >- farm2;
                farm2.trading_partners >- farm1;
            }
            do die;
        }
        
        // Phase 2: Edge Formation
        // Use targeted sampling to maintain computational efficiency
        list<farm> active_farms <- farm where (each.status != "removed");
        int n_active <- length(active_farms);
        
        if n_active > 1 {
            // Calculate current mean degree
            list<int> current_degrees <- active_farms collect length(each.trading_partners);
            current_mean_degree <- mean(current_degrees);
            
            // Calculate how many edges to form based on mean degree target
            int current_edges <- length(edge_connection);
            float desired_edges <- (n_active * target_mean_degree) / 2.0;
            int edges_to_form <- int(desired_edges - current_edges);
            
            // Apply mean degree correction - more aggressive for sparse networks
            float degree_ratio <- target_mean_degree > 0 ? current_mean_degree / target_mean_degree : 0;
            
            if degree_ratio > 1.2 {
                edges_to_form <- int(edges_to_form * 0.7); // Slightly slow down if above target
            } else if degree_ratio < 0.9 {
                edges_to_form <- int(edges_to_form * 2.0); // Aggressively speed up if below target
            }
            
            // Ensure minimum edge formation when below target
            if current_mean_degree < target_mean_degree * 0.9 {
                edges_to_form <- max(edges_to_form, int(n_active * 0.1)); // At least 10% of nodes get new edges
            }
            
            if edges_to_form > 0 {
                // Sample potential edges to form
                int attempts <- 0;
                int max_attempts <- min(edges_to_form * 10, n_active * 5);
                int formed <- 0;
                
                loop while: formed < edges_to_form and attempts < max_attempts {
                    attempts <- attempts + 1;
                    
                    // Random selection of farms
                    farm f1 <- one_of(active_farms);
                    farm f2 <- one_of(active_farms);
                    
                    if f1 = f2 or f1.trading_partners contains f2 {
                        continue;
                    }
                    
                    // Calculate formation probability
                    float log_odds <- calculate_formation_log_odds(f1, f2, false);
                    float prob_form <- 1.0 / (1.0 + exp(-log_odds));
                    
                    // Apply mean degree correction - more aggressive
                    float degree_correction <- 1.0;
                    if current_mean_degree > target_mean_degree * 1.2 {
                        degree_correction <- 0.7;
                    } else if current_mean_degree < target_mean_degree * 0.9 {
                        degree_correction <- 2.0;  // Double the probability when below target
                    }
                    prob_form <- min(0.98, prob_form * degree_correction);  // Increased max from 0.95
                    
                    if flip(prob_form) {
                        create edge_connection {
                            farm1 <- f1;
                            farm2 <- f2;
                            formation_time <- current_week;
                            edge_duration <- 0;
                        }
                        
                        f1.trading_partners << f2;
                        f2.trading_partners << f1;
                        formed <- formed + 1;
                        current_edges <- current_edges + 1;
                    }
                }
            }
        }
        
        // Update network metrics
        total_edges <- length(edge_connection);
        if n_active > 1 {
            current_density <- (2.0 * total_edges) / (n_active * (n_active - 1));
            list<int> final_degrees <- active_farms collect length(each.trading_partners);
            current_mean_degree <- mean(final_degrees);
        } else {
            current_density <- 0.0;
            current_mean_degree <- 0.0;
        }
    }
    
    action compute_network_statistics {
        // Compute various network statistics for STERGM monitoring
        list<farm> active_farms <- farm where (each.status != "removed");
        int n <- length(active_farms);
        
        if n = 0 {
            network_statistics <- map([]);
            return;
        }
        
        // Degree distribution statistics
        list<int> degrees <- active_farms collect length(each.trading_partners);
        float mean_deg <- mean(degrees);
        float variance_deg <- variance(degrees);
        int max_deg <- max(degrees);
        
        // Homophily measure (proportion of same-type edges)
        float homophily <- 0.0;
        int same_type_edges <- 0;
        int total_counted <- 0;
        
        ask edge_connection {
            if farm1.status != "removed" and farm2.status != "removed" {
                total_counted <- total_counted + 1;
                if farm1.farm_type = farm2.farm_type {
                    same_type_edges <- same_type_edges + 1;
                }
            }
        }
        
        if total_counted > 0 {
            homophily <- same_type_edges / total_counted;
        }
        
        // Clustering coefficient (transitivity)
        float clustering <- calculate_clustering_coefficient(active_farms);
        
        // Store statistics
        network_statistics <- [
            "nodes"::float(n),
            "edges"::float(total_edges),
            "density"::current_density,
            "mean_degree"::mean_deg,
            "variance_degree"::variance_deg,
            "max_degree"::float(max_deg),
            "homophily"::homophily,
            "clustering"::clustering,
            "isolates"::float(length(active_farms where (length(each.trading_partners) = 0)))
        ];
        
        // Track history
        density_history << current_density;
        transitivity_history << clustering;
    }
    
    float calculate_clustering_coefficient(list<farm> farms) {
        // Calculate global clustering coefficient
        float total_triangles <- 0.0;
        float total_triplets <- 0.0;
        
        ask farms {
            int k <- length(trading_partners);
            if k >= 2 {
                // Count triangles
                int triangles <- 0;
                loop i from: 0 to: k - 2 {
                    loop j from: i + 1 to: k - 1 {
                        farm neighbor1 <- trading_partners[i];
                        farm neighbor2 <- trading_partners[j];
                        if neighbor1.trading_partners contains neighbor2 {
                            triangles <- triangles + 1;
                        }
                    }
                }
                total_triangles <- total_triangles + triangles;
                total_triplets <- total_triplets + (k * (k - 1) / 2.0);
            }
        }
        
        if total_triplets > 0 {
            return total_triangles / total_triplets;
        } else {
            return 0.0;
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
        ask newly_infected {
            status <- "infected";
            infection_time <- current_week;
            infected_by_type[farm_type] <- infected_by_type[farm_type] + 1;
        }
        
        num_infected_farms <- num_infected_farms + length(newly_infected);
        num_susceptible_farms <- num_susceptible_farms - length(newly_infected);
        
        if first_infection_week < 0 {
            first_infection_week <- current_week;
        }
    }
    
    action select_index_case {
        // Select index case from medium farms with high connectivity
        list<farm> medium_farms <- farm where (each.farm_type = 2);
        
        if !empty(medium_farms) {
            // Prefer well-connected farms
            medium_farms <- medium_farms sort_by (-length(each.trading_partners));
            farm index_farm <- first(medium_farms);
            
            if length(index_farm.trading_partners) = 0 {
                // If no connected medium farm, choose randomly
                index_farm <- one_of(medium_farms);
            }
            
            index_farm.status <- "infected";
            index_farm.infection_time <- 0;
            
            num_infected_farms <- 1;
            num_susceptible_farms <- num_susceptible_farms - 1;
            infected_by_type[2] <- 1;
            first_infection_week <- 0;
            
            index_id <- index_farm.farm_id;
            
            write "Index case: " + index_id + " with " + 
                  length(index_farm.trading_partners) + " connections";
        }
    }
    
    action apply_culling {
        list<farm> infected_list <- farm where (each.status = "infected");
        int num_culled <- length(infected_list);
        
        if num_culled = 0 { return; }
        
        write "=== CULLING INTERVENTION at Week " + current_week + " ===";
        write "Culling " + num_culled + " infected farms";
        
        ask infected_list {
            culled_by_type[farm_type] <- culled_by_type[farm_type] + 1;
            infected_by_type[farm_type] <- max(0, infected_by_type[farm_type] - 1);
            status <- "removed";
        }
        
        // Remove all edges connected to culled farms
        list<edge_connection> edges_to_remove <- edge_connection where 
            (infected_list contains each.farm1 or infected_list contains each.farm2);
        
        ask edges_to_remove { do die; }
        
        ask infected_list {
            ask trading_partners { trading_partners >- myself; }
            trading_partners <- [];
        }
        
        num_removed_farms <- num_removed_farms + num_culled;
        num_infected_farms <- 0;
    }
    
    action output_initial_state {
        write "========================================";
        write "=== STERGM MODEL INITIALIZATION ===";
        write "========================================";
        write "Total farms: " + farms_count;
        write "";
        write "PRIMARY FIXED PARAMETERS:";
        write "  - TARGET MEAN DEGREE: " + target_mean_degree + 
              " (" + connectivity_scenario + " scenario)";
        write "    * Low: 0.5, Baseline: 0.75, High: 1.0";
        write "  - Transmission probability: " + transmission_probability;
        write "  - Contact rates: Fixed matrix";
        write "  - Edge duration: " + (edge_duration_enabled ? 
              "ENABLED (min " + edge_duration + " weeks)" : 
              "DISABLED (departure rate: " + departure_rate + ")");
        write "";
        write "Network Statistics:";
        write "  - Actual mean degree: " + round(current_mean_degree * 100) / 100;
        write "  - Target density: " + round(target_density * 10000) / 10000;
        write "  - Actual density: " + round(current_density * 10000) / 10000;
        write "  - Total edges: " + total_edges + " (target: " + target_edges + ")";
        write "  - Clustering: " + round(network_statistics["clustering"] * 1000) / 1000;
        write "  - Homophily: " + round(network_statistics["homophily"] * 100) / 100 + "%";
        write "";
        write "STERGM Parameters:";
        write "  - Formation intercept: " + theta_formation_edges;
        write "  - Formation homophily: " + theta_formation_homophily;
        write "  - Dissolution intercept: " + theta_dissolution_edges;
        write "  - Dissolution duration: " + theta_dissolution_duration;
        write "";
        write "Intervention:";
        write "  - Culling: " + (culling_enabled ? "ON at week " + culling_timing : "OFF");
        write "========================================";
    }
    
    action finalize {
        do compute_network_statistics;
        
        int total_affected <- num_removed_farms + num_infected_farms;
        float attack_rate <- (total_affected * 100.0) / farms_count;
        
        write "";
        write "========================================";
        write "=== STERGM SIMULATION COMPLETE ===";
        write "========================================";
        write "Duration: " + current_week + " weeks";
        write "";
        write "FINAL EPIDEMIC STATUS:";
        write "  Susceptible: " + num_susceptible_farms;
        write "  Infected: " + num_infected_farms;
        write "  Removed: " + num_removed_farms;
        write "  Attack rate: " + round(attack_rate) + "%";
        write "";
        write "FINAL NETWORK STATISTICS:";
        write "  Density: " + round(current_density * 1000) / 1000;
        write "  Mean degree: " + round(network_statistics["mean_degree"] * 100) / 100;
        write "  Clustering: " + round(network_statistics["clustering"] * 1000) / 1000;
        write "  Homophily: " + round(network_statistics["homophily"] * 100) + "%";
        write "";
        write "NETWORK EVOLUTION:";
        write "  Mean density: " + round(mean(density_history) * 1000) / 1000;
        write "  Density variance: " + round(variance(density_history) * 10000) / 10000;
        write "========================================";
        
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
                          target_mean_degree + "," + 
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