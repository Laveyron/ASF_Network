model ASF_Network_Simulation

global {
    action output_initial_state {
        do calculate_stats();
        
        write "========================================";
        write "=== INITIAL STATE SUMMARY ===";
        write "========================================";
        write "Total farms: " + length(farm) + " (Small:" + farms_by_type[1] + 
              " Medium:" + farms_by_type[2] + " Large:" + farms_by_type[3] + ")";
        write "Network: " + total_edges + " edges, mean degree: " + 
              overall_mean_degree;
//        write "Isolated farms: " + isolated + " (" + round(isolated * 100.0 / farms_count) + "%)";
        write "Index case selected: Medium farm " + index_id;
        write "Parameters: transmission probability = " + transmission_probability;
        write "Intervention: culling " + (culling_enabled ? "ON at week " + culling_timing : "OFF");
//        write total_edges;
//        write mean_degree;
        if total_edges < (target_edges * 0.8) {
            write "WARNING: Network formation achieved only " + 
                  round(total_edges * 100.0 / target_edges) + "% of target connectivity";
        }
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
}