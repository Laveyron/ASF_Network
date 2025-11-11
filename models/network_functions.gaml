model ASF_Network_Simulation

global {
	action set_connectivity_scenario {
        float scenario_value;
        if connectivity_scenario = "Thấp" {
            scenario_value <- 0.5;
        } else if connectivity_scenario = "Cao" {
            scenario_value <- 1.0;
        } else {
            scenario_value <- 0.75;
        }
        
        overall_mean_degree <- scenario_value;
        target_edges <- int((length(farm) * overall_mean_degree) / 2.0);
    }
    
    float calculate_distance_kernel(farm f1, farm f2) {
        float dist_km <- f1 distance_to f2 / 1000.0;
        float weight <- k0 / (1.0 + (dist_km / r0) ^ alpha);
        return weight;
    }
    
    action initialize_network {
        list<farm> all_farms <- list(farm);
        farms_count <- length(all_farms);
        
        int edges_created <- 0;
        
        loop while: (edges_created < target_edges) {     
            farm f1 <- one_of(all_farms);
            farm f2 <- one_of(all_farms);
            
            if (f1.farm_id = f2.farm_id) or (f1.trading_partners contains f2) { continue; }
            
            if length(f1.trading_partners) >= 5 or length(f2.trading_partners) >= 5 {
                continue;
            }
            
            string key <- string(f1.farm_type) + "-" + string(f2.farm_type);
            float base_prob <- contact_rates contains_key key ? contact_rates[key] : 0.0;
            
            if base_prob <= 0 { continue; }
            
            float distance_weight <- calculate_distance_kernel(f1, f2);
            float final_prob <- base_prob * distance_weight;
            
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
        }
        
        total_edges <- edges_created;
        
        if !empty(all_farms) {
            float actual_mean_degree <- mean(all_farms collect float(length(each.trading_partners)));
            mean_degree <- actual_mean_degree;
        }
        
        isolated <- length(all_farms where (length(each.trading_partners) = 0));
    }
    
    action update_network_stergm {
        // Dissolution Phase
        list<edge_connection> edges_to_remove;

        edges_to_remove <- edge_connection where 
            (each.farm1 = nil or each.farm2 = nil or 
             each.farm1.status = "removed" or each.farm2.status = "removed");

        if edge_duration_enabled {
            edges_to_remove <- edges_to_remove + (edge_connection where 
                (each.farm1.status != "removed" and each.farm2.status != "removed" and
                 (current_week - each.formation_time) >= each.expected_duration));
        } else {
            edges_to_remove <- edges_to_remove + (edge_connection where 
                (each.farm1.status != "removed" and each.farm2.status != "removed" and
                 flip(dissolution_rate)));
        }

        if !empty(edges_to_remove) {
            ask edges_to_remove {
                if farm1 != nil and farm2 != nil {
                    farm1.trading_partners >- farm2;
                    farm2.trading_partners >- farm1;
                }
                do die;
            }
        }

        // Formation Phase
        int current_edges <- length(edge_connection);
        int edges_needed <- max(0, target_edges - current_edges);

        if edges_needed > 0 {
            int formed <- form_edges(edges_needed);
        }
    }   

    int form_edges(int target_new_edges) {
        list<farm> active_farms <- farm where (each.status != "removed");
        if length(active_farms) < 2 { return 0; }
        
        int formed <- 0;
        
        int max_degree <- 5;
        list<farm> available_farms <- active_farms where (length(each.trading_partners) < max_degree);
        
        if length(available_farms) < 2 { 
            return 0;
        }
        
        loop while: formed < target_new_edges {
            farm f1 <- one_of(available_farms);
            farm f2 <- one_of(available_farms);
            
            if f1 = f2 { continue; }
            if f1.trading_partners contains f2 { continue; }
            
    
            if length(f1.trading_partners) >= max_degree or length(f2.trading_partners) >= max_degree {
                available_farms >- f1;
                available_farms >- f2;
                if length(available_farms) < 2 { break; }
                continue;
            }
            
            string key <- string(f1.farm_type) + "-" + string(f2.farm_type);
            float base_prob <- contact_rates contains_key key ? contact_rates[key] : 0.0;
            
            if base_prob <= 0 { continue; }
            
            float distance_weight <- calculate_distance_kernel(f1, f2);
            float stergm_prob <- base_prob * distance_weight;
            
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
                
                if length(f1.trading_partners) >= max_degree {
                    available_farms >- f1;
                }
                if length(f2.trading_partners) >= max_degree {
                    available_farms >- f2;
                }
                
                if length(available_farms) < 2 { break; }
            }
        }
        
        return formed;
    }
}