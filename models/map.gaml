model ASF_Network_Simulation

species boundary {
    aspect default {
        draw shape color: rgb(207, 207, 207) border: rgb(37, 37, 37) width: 1;
    }
}

species road {
    aspect default {
        draw shape color: rgb(60, 60, 60) width: 0.5;
    }
}

species water {
    aspect default {
        draw shape color: rgb(120, 172, 206) border: rgb(120, 172, 206);
    }
}

species waterways {
    aspect default {
        draw shape color: rgb(99, 140, 201) border: rgb(99, 140, 201);
    }
}

global {
    action load_spatial_data {
        create boundary from: boundary_file;
        create road from: roads_file;
        create water from: water_file;
        create waterways from: waterways_file;
    }
}