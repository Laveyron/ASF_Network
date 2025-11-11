model ASF_Network_Simulation

global {
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
}