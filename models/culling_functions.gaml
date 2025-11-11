model ASF_Network_Simulation

global {
    action apply_culling {
        list<farm> infected_list <- farm where (each.status = "infected");
        int num_culled <- length(infected_list);
        
        if num_culled = 0 { return; }
        
        write "";
        write "========================================";
        write "=== CULLING INTERVENTION ===";
        write "========================================";
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
}