model ASF_Network_Simulation

global {
    action calculate_stats {
        list<farm> active <- farm where (each.status != "removed");
        
        if !empty(active) {
            list<int> degrees <- active collect length(each.trading_partners);
            mean_degree <- mean(degrees);
        } else {
            mean_degree <- 0.0;
        }
        
        total_edges <- length(edge_connection);
        
        isolated <- length(farm where (length(each.trading_partners) = 0));
    }
    
    action update_statistics {
        num_infected_farms <- length(farm where (each.status = "infected"));
        num_susceptible_farms <- length(farm where (each.status = "susceptible"));
        num_removed_farms <- length(farm where (each.status = "removed"));
        
        num_small_infected <- infected_by_type[1];
        num_medium_infected <- infected_by_type[2];
        num_large_infected <- infected_by_type[3];
        
        num_small_culled <- culled_by_type[1];
        num_medium_culled <- culled_by_type[2];
        num_large_culled <- culled_by_type[3];
    }
}