function [mod, top_of_specimen] = fn_set_ply_material(mod, matl_1, matl_2, specimen_brdy_pts, n_layers, n_plys_per_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Calculate ply layer heights so ply material is applied to specific layers
sbp = specimen_brdy_pts; %tmp for readability
specimen_height = sbp(3,2) - sbp(1,2);
ply_height = specimen_height/32;
el_height = mod.el_height; %Calculated in fn_isometric_structured mesh
height_offset = rem(ply_height,el_height);
%Alternate between upper and lower ply heights
ply_height_lower = ply_height - height_offset;
ply_height_upper = ply_height + (el_height - height_offset);

%Calculate specimen height so that it is a multiple of element layers
specimen_offset = rem(sbp(1,2),el_height);
specimen_height_from_bottom = sbp(1,2) - specimen_offset;
specimen_width = specimen_brdy_pts(2,1);
%Apply margin to x so that -ve x values are captured
safety_margin_perc = 0.1;
safety_margin_x = specimen_width * safety_margin_perc;

%Assign materials
matl_oscillator = 1;
height_completed = specimen_height_from_bottom;
if rem(n_layers,2)
    error("n_layers must be even")
end
%Loop over ply types
for ply_type = 1:n_layers/n_plys_per_type
    %This changes ply material when needed
    if matl_oscillator > 0
        %first material
        matl_i = matl_1;
    else
        %second material
        matl_i = matl_2;
    end
    %Loop over layers in stack
    for layer_in_type = 1:n_plys_per_type
        target_layer = (ply_type - 1)*n_plys_per_type + layer_in_type;
        %Stay close to true ply layer height using upper and lower heights
        true_height = specimen_height_from_bottom + (target_layer - 1)*ply_height;
        if height_completed > true_height
            ply_height_target = ply_height_lower;
        else
            ply_height_target = ply_height_upper;
        end
                
        %Boundary points of target layer
        bdry_pts = [0,               0
                    specimen_width,  0
                    specimen_width,  ply_height_target
                    0,               ply_height_target];
        %Translate upwards by height completed (starts from bottom of specimen)
        bdry_pts = bdry_pts + height_completed*[0 1
                                                0 1
                                                0 1
                                                0 1];
        %Apply safety margin
        bdry_pts = bdry_pts + safety_margin_x*[-1 0
                                                1 0
                                                1 0
                                               -1 0];
        %Change materials
        mod = fn_set_els_inside_bdry_to_mat(mod, bdry_pts, matl_i);
        %Update height of material changed in specimen
        height_completed = height_completed + ply_height_target;
    end
    %Flip material type after n_plys_per_type loops completed
    matl_oscillator = matl_oscillator*-1;
    %CFRP specimen is symmetric
    if target_layer == n_layers/2
        matl_oscillator = matl_oscillator*-1;
    end
end

%Return top of specimen (because new height != defined specimen height)
top_of_specimen = height_completed;

end