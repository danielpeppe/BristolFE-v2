function [mod, top_of_specimen] = fn_set_ply_material(mod, op, matl_1, matl_2, specimen_brdy_pts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Extract options
n_ply_layers = op.n_ply_layers;
n_plys_per_type = op.n_plys_per_type;
ply_symmetry = op.ply_symmetry;
upper_water_present = op.upper_water_present;
ply_type_boundary = op.ply_type_boundary;

%Calculate ply layer heights so ply material is applied to specific layers
sbp = specimen_brdy_pts; %tmp for readability
specimen_height = sbp(3,2) - sbp(1,2);
ply_height = specimen_height/n_ply_layers;
el_height = mod.el_height; %Calculated in fn_isometric_structured mesh
height_offset = rem(ply_height,el_height);
%Alternate between upper and lower ply heights
ply_height_lower = ply_height - height_offset;
ply_height_upper = ply_height + (el_height - height_offset);

%Calculate specimen height so that it is a multiple of element layers
specimen_offset = rem(sbp(1,2),el_height);
specimen_height_from_bottom = sbp(1,2) - specimen_offset;
specimen_width = specimen_brdy_pts(2,1);
%Apply margin to x so that -ve x values are captured, also ensures last
%layer leaves no water if upper water layer is present
safety_margin_perc = 0.3;
safety_margin = specimen_width * safety_margin_perc;

%Assign materials
matl_oscillator = 1;
height_completed = specimen_height_from_bottom;
if rem(n_ply_layers,2)
    error("n_ply_layers must be even")
end

%Loop over ply types
for ply_type = 1:n_ply_layers/n_plys_per_type
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
            ply_target_height = ply_height_lower;
        else
            ply_target_height = ply_height_upper;
        end
        %Check if target_layer is final layer if upper water is present (INEFFICIENT CODE, BUT DOESNT MATTER)
        if (target_layer == n_ply_layers) && ~upper_water_present
            ply_target_height = (1 + safety_margin_perc)*ply_target_height;
        end
        %Set materials to target layer
        [mod, height_completed] = set_target_layer_material(mod, matl_i, specimen_width, ply_target_height, height_completed, safety_margin, 0);
        %Update height of material changed in specimen
        height_completed = height_completed + ply_target_height;
        %Set inter-ply-type inertia boundary for reflections
        if ply_type_boundary && (layer_in_type < n_plys_per_type || (ply_symmetry && (target_layer == n_ply_layers/2)))
            mod = set_target_layer_material(mod, matl_i + 1, specimen_width, mod.el_height, height_completed, safety_margin, 1);
        end

    end
    %Flip material type after n_plys_per_type loops completed
    matl_oscillator = matl_oscillator*-1;
    %CFRP specimen is symmetric
    if ply_symmetry && (target_layer == n_ply_layers/2)
        matl_oscillator = matl_oscillator*-1;
    end
end

%Calculate deviation from true specimen height
model_specimen_height = height_completed - sbp(1,2);
model_accuracy = round(100*model_specimen_height/specimen_height, 2);
fprintf("Specimen model accuracy: %.2f%%\n", model_accuracy)
if abs(model_accuracy) < 95
    disp('CAUTION: Specimen model accuracy is less than 95%')
end
%Return top of specimen (because new height != defined specimen height)
top_of_specimen = height_completed;
end


function [mod, height_completed] = set_target_layer_material(mod, matl_i, specimen_width, target_height, height_completed, safety_margin, inter_ply_boundary)
%Boundary points of target layer
bdry_pts = [0,               0
            specimen_width,  0
            specimen_width,  target_height
            0,               target_height];
%Translate upwards by height completed (starts from bottom of specimen)
bdry_pts = bdry_pts + height_completed*[0 1
                                        0 1
                                        0 1
                                        0 1];
%Apply safety margin
bdry_pts = bdry_pts + safety_margin*[-1 0
                                      1 0
                                      1 0
                                     -1 0];

if inter_ply_boundary
    if target_height ~= mod.el_height
        error('Ply boundary material must be 1 layer thick')
    end
    bdry_pts = bdry_pts - mod.el_height*[0 1
                                         0 1
                                         0 1
                                         0 1];
end
%Change materials
mod = fn_set_els_inside_bdry_to_mat(mod, bdry_pts, matl_i);
end