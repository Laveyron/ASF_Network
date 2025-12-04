model ASF_Network_Simulation

global {
    list<farm> do_transmission {
    	list<farm> pending_infection <- [];
    	list<farm> infected_list <- farm where (each.status = "infected");
    	float weight <- 1.0;

    	ask infected_list {
        loop partner over: trading_partners {
            if partner.status = "susceptible" {
                string key <- string(self.farm_type) + "-" + string(partner.farm_type);
                float contact_rate <- contact_rates contains_key key ? contact_rates[key] : 0.0;
                
                if spatial_sampling_enabled {
                	float dist_km <- (self distance_to partner) / 1000.0;
        			weight <- k0 / (1.0 + ((dist_km / r0) ^ alpha));
                }
                
                if contact_rate > 0 {
                	int actual_contacts <- poisson(contact_rate);
                	
                	if (actual_contacts > 0) {
                		loop i from: 1 to: actual_contacts {
                			
//                			write transmission_probability*weight;
                			
                			if flip(transmission_probability * weight) {
                				pending_infection << partner;
                				break;
                			}
                		} 
                	}
                }
            }
        }
    }

    
        return remove_duplicates(pending_infection);
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
}