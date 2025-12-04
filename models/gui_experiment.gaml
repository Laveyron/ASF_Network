model ASF_Network_Simulation

import "main.gaml"

experiment ASF_GUI_Simulation type: gui {
    parameter "Mức độ kết nối" var: connectivity_scenario among: ["Thấp", "Cơ sở", "Cao"];
    parameter "Xác suất truyền nhiễm" var: transmission_probability min: 0.1 max: 1.0 step: 0.1;
    parameter "Kích hoạt biện pháp tiêu hủy" var: culling_enabled;
    parameter "Tuần tiêu hủy" var: culling_timing min: 1 max: 52 step: 1;
    parameter "Kích hoạt duy trì liên kết" var: edge_duration_enabled;
    parameter "Thời gian duy trì liên kết (tuần)" var: edge_duration min: 6 max: 52 step: 1;
    parameter "Áp dụng hàm kernel truyền nhiễm" var: spatial_sampling_enabled init: false;
    parameter "Experiment ID" var: experiment_id <- "";
    
    reflex capture when: mod(cycle, 1) = 0 {
        ask simulations {
            save (snapshot(self, "map", {1920.0, 1080.0})) to: "../includes/output/map/" + experiment_id + "map" + string(cycle) + ".png";
            save (snapshot(self, "Infections by Farm Type (Pie)", {800.0, 600.0})) to: "../includes/output/map/"+ experiment_id + "pie" + string(cycle) + ".png";
            save (snapshot(self, "Infection Trends Over Time", {1200.0, 600.0})) to: "../includes/output/map/"+ experiment_id + "line" + string(cycle) + ".png";
        }
    }
    
    output synchronized: true {
        display map type: opengl 
        background: rgb(245, 247, 250)
        antialias: true
        axes: false {

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
            font header_reg <- font("Arial", fh - 2, #bold);
            font body_bold <- font("Arial", fb, #bold);
            font body_reg <- font("Arial", fb);
            font small_bold <- font("Arial", fs, #bold);
            font small_reg <- font("Arial", fs);

            // Layout metrics
            float sh <- h * 0.1;
            float ls <- h * 0.06;
            float y <- my;

            // Header
            draw "ASF STATUS" at: {cx, y + sh/2}
                color: #black font: header_reg anchor: #center;
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
            draw string(length(farm)) at: {rx - mx + 500, y + fh2/2}
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
        
        display "Infections by Farm Type (Pie)" type: java2D {
            chart "Tỷ lệ nhiễm bệnh theo loại hình trang trại" 
                type: pie 
                background: #white
                color: #black {
                    data "Trang trại nhỏ" value: num_small_infected color: rgb(102, 194, 165);
                    data "Trang trại vừa" value: num_medium_infected color: rgb(252, 141, 98);
                    data "Trang trại lớn" value: num_large_infected color: rgb(141, 160, 203);
                }
        }
        
        display "Infection Trends Over Time" type: java2D {
            chart "Diễn biến số trang trại nhiễm bệnh theo loại hình trang trại" 
                type: series
                style: stack
                background: #white
                color: #black
                x_label: "Tuần"
                y_label: "Số trang trại nhiễm bệnh"
                x_range: [0, max_simulation_weeks]
                x_tick_unit: 5 {
                    data "Trang trại nhỏ" value: num_small_infected color: rgb(102, 194, 165) thickness: 1;
                    data "Trang trại vừa" value: num_medium_infected color: rgb(252, 141, 98) thickness: 1;
                    data "Trang trại lớn" value: num_large_infected color: rgb(141, 160, 203) thickness: 1;
                    data "Tổng" value: num_infected_farms color: rgb(50, 50, 50) marker: false thickness: 2;
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
