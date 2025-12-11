model ASF_Network_Simulation

global {
    // ERGM Parameters
    float theta_edges <- -3.0;   
    float theta_nodematch <- 1.4;
    float theta_degree <- -0.05;
    float theta_isolates <- -5.5;
    
    map<farm, list<farm>> potential_partners <- [];
    map<farm, list<farm>> neighbors <- [];
    bool spatial_sampling_enabled <- false;
    int proposals;
    
    action set_connectivity_scenario {
        float scenario_value;
        if connectivity_scenario = "Thấp" {
            scenario_value <- 0.5;
            proposals <- 4480;
        } else if connectivity_scenario = "Cao" {
            scenario_value <- 1.0;
            proposals <- 10280;
        } else {
            scenario_value <- 0.75;
            proposals <- 7100;
        }
        
        overall_mean_degree <- scenario_value;
        target_edges <- int((length(farm) * overall_mean_degree) / 2.0);
    }
    
    // ERGM Log-ratio
    float compute_formation_log_ratio(farm f1, farm f2) {
        int deg1 <- length(f1.trading_partners);
        int deg2 <- length(f2.trading_partners);
        
        // Change statistics
        int delta_nodematch <- (f1.farm_type = f2.farm_type) ? 1 : 0;
        int delta_isolates <- ((deg1 = 0) ? -1 : 0) + ((deg2 = 0) ? -1 : 0);
        int delta_degree <- 2;
        int delta_edges <- 1;
        
        // ERGM log-odds
        float log_odds <- theta_edges * delta_edges + 
                         theta_nodematch * delta_nodematch + 
                         theta_degree * delta_degree + 
                         theta_isolates * delta_isolates;
        
        return log_odds;
    }
    
    // Metropolis-Hastings acceptence
    float mh_acceptance_probability(float log_ratio) {
        if log_ratio <= -10.0 { return 0.0; }
        if log_ratio >= 5.0 { return 1.0; }
        return min(1.0, exp(log_ratio));
    }
    
    // Dyad sampling
    list<farm> sample_dyad(list<farm> active) {
        farm f1 <- one_of(active);
        farm f2 <- one_of(potential_partners[f1]);
        return [f1, f2];
    }
    
    action init_neighbors {
        ask farm {
            neighbors[self] <- farm where (
                each != self and 
                ((each distance_to self) / 1000.0) <= 3.0
            );
        }
    }
    
    // Dyad spatial sampling
    list<farm> spatial_dyad(list<farm> active) {
        farm f1 <- one_of(active);
        farm f2;
        
        list<farm> nearby <- neighbors[f1] where (each.status != "removed");
    	float spatial_probability <- 0.95;
    
    	if !empty(nearby) and flip(spatial_probability) {
        	f2 <- one_of(nearby);
    	} else {
        	f2 <- one_of(active where (each != f1));
    	}
            
        return [f1, f2];
    }
    
    // ERGM network initialization
    action initialize_network {
        list<farm> all_farms <- list(farm);
        farms_count <- length(all_farms);
        
        if (spatial_sampling_enabled) {
        	do init_neighbors();
        }
        
        potential_partners <- [];
        ask all_farms {
            potential_partners[self] <- all_farms - [self];
        }
        
        
        loop step from: 1 to: proposals {
        	list<farm> dyad;
    
    		if (spatial_sampling_enabled) {
        		dyad <- spatial_dyad(all_farms);
    		} else {
        		dyad <- sample_dyad(all_farms);
    		}
            
        	farm f1 <- dyad[0];
        	farm f2 <- dyad[1];
            
            if f1 = nil or f2 = nil { continue; }
            
            if f1.trading_partners contains f2 { continue; }
            
            float log_ratio <- compute_formation_log_ratio(f1, f2);
                
            if flip(mh_acceptance_probability(log_ratio)) {
            	create edge_connection {
                	farm1 <- f1;
                    farm2 <- f2;
                  	formation_time <- 0;
                    expected_duration <- edge_duration;
                }
                f1.trading_partners << f2;
                f2.trading_partners << f1;
            }
        }
    }
    
    // STERGMs
    action update_network_stergm {
        
        // Dissolution
        list<edge_connection> edges_to_remove <- [];
        
        ask edge_connection {
            bool remove <- false;
            
            if farm1 = nil or farm2 = nil or 
               farm1.status = "removed" or farm2.status = "removed" {
                remove <- true;
            }
            else if edge_duration_enabled {
                remove <- (current_week - formation_time) >= expected_duration;
            }
            else {
                remove <- flip(dissolution_rate);
            }
            
            if remove { edges_to_remove << self; }
        }
        
        int dissolved <- length(edges_to_remove);
        
        ask edges_to_remove {
            if farm1 != nil and farm2 != nil {
                farm1.trading_partners >- farm2;
                farm2.trading_partners >- farm1;
            }
            do die;
        }
        
        // Formation
        list<farm> active <- farm where (each.status != "removed");
        int n_active <- length(active);
        
        if n_active < 2 { return; }
        
        if edge_duration_enabled {
        	ask edge_connection {
				if (current_week - formation_time) < expected_duration {
        			return;
        		}
        	}
        }
        
        if empty(potential_partners) {
            potential_partners <- [];
            ask active {
                potential_partners[self] <- active - [self];
            }
        }
        
        loop i from: 1 to: proposals {
            list<farm> dyad;
    
    		if (spatial_sampling_enabled) {
        		dyad <- spatial_dyad(active);
    		} else {
        		dyad <- sample_dyad(active);
    		}
            
            farm f1 <- dyad[0];
            farm f2 <- dyad[1];
            
            if f1 = nil or f2 = nil { continue; }
            
            if f1.trading_partners contains f2 { continue; }
            

            float log_ratio <- compute_formation_log_ratio(f1, f2);
            
            if flip(mh_acceptance_probability(log_ratio)) {
                create edge_connection {
                    farm1 <- f1;
                    farm2 <- f2;
                    formation_time <- current_week;
                    expected_duration <- edge_duration;
                }
                f1.trading_partners << f2;
                f2.trading_partners << f1;
            }
        }
    }
}