function mod = fn_add_ply_waviness(mod, op, comp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Create ply layer probability distribution
specimen_width = comp.specimen_width;
safety_margin_x = specimen_width * 0.3; %arbitrary (just to capture all x-nodes)


for i = 1:op.n_ply_layers
    upper_layer_height = comp.ply_location_tracker{i, 1};
    lower_layer_height = comp.ply_location_tracker{i, 2};
    
    upper_loc1 = [upper_layer_height (0 - safety_margin_x)];
    upper_loc2 = [upper_layer_height (specimen_width + safety_margin_x)];
    lower_loc1 = [lower_layer_height (0 - safety_margin_x)];
    lower_loc2 = [lower_layer_height (specimen_width + safety_margin_x)];
    
    %[y x]
    % el_size = mod.el_height / sqrt(3);
    nds_upper = fn_find_nodes_on_line(mod.nds, upper_loc1, upper_loc2, el_size / 2);
    nds_lower = fn_find_nodes_on_line(mod.nds, lower_loc1, lower_loc2, el_size / 2);

    mod.nds(nodes_upper, 2) = mod.nds(nds_upper, 2);

end

end

