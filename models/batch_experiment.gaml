model ASF_Network_Simulation

experiment Batch_Connectivity_Scenarios type: batch repeat: 20 keep_seed: false parallel: 20 until: (current_week >= max_simulation_weeks) {
    
    parameter "Enable Culling" var: culling_enabled init: true;
    parameter "Culling Week" var: culling_timing init: 6;
    parameter "Enable Edge Duration" var: edge_duration_enabled init: false;
    parameter "Edge Duration" var: edge_duration init: 26;
    parameter "Transmission Probability" var: transmission_probability init: 0.6;
    parameter "Connectivity" var: connectivity_scenario init: "Cơ sở";
}