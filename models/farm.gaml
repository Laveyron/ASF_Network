model ASF_Network_Simulation

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
            match 3 { size <- 170.0; }
        }
        
        rgb color_farm;
        switch status {
            match "susceptible" { color_farm <- rgb(27, 158, 119); }
            match "infected" { color_farm <- rgb(255, 0, 255); size <- 120.0;}
            match "removed" { color_farm <- rgb(255, 0, 0); }
        }
        
        draw circle(size) color: color_farm border: #black width: 0.1;
        
//        if status = "infected" {
//        	draw circle (3000) color: rgb(255, 0, 0, 50) border: #red;
//        }
    }
    
    aspect network_view {
        float size;
        switch farm_type {
            match 1 { size <- 20.0; }
            match 2 { size <- 60.0; }
            match 3 { size <- 100.0; }
        }
        
        rgb color_farm;
        switch status {
            match "susceptible" { color_farm <- rgb(27, 158, 119); }
            match "infected" { color_farm <- rgb(255, 0, 255); size <- 200.0;}
            match "removed" { color_farm <- rgb(255, 0, 0); }
        }
        
        draw circle(size) color: color_farm border: #black;
        
//        draw circle (550) color: rgb(255, 0, 0, 50) border: #red;
       
        loop partner over: trading_partners {
            draw line([location, partner.location]) color: #black width: 1;
        }
    }
}

global {
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
        }
    }
    
}