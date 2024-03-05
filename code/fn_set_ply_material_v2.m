function [mod, top_of_specimen] = fn_set_ply_material_v2(mod, op, matls, specimen_brdy_pts, model_height)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fprintf("Set Ply Materials (v2)\n")
%Ensure specimen_brdy_pts is a double precision matrix
specimen_brdy_pts = double(specimen_brdy_pts);

%Put key geometry consts in a struct
geom.model_height = model_height;
geom.specimen_width = specimen_brdy_pts(2,1);
safety_margin_perc = 0.3; %arbitrary
geom.safety_margin = geom.specimen_width * safety_margin_perc; %Apply margin to x so that -ve x values are captured
geom.specimen_height = specimen_brdy_pts(3,2) - specimen_brdy_pts(1,2);

%Define and declare options used
n_ply_layers = op.n_ply_layers;
n_plys_per_type = op.n_plys_per_type;
ply_symmetry = op.ply_symmetry;
interply_boundary = op.interply_boundary;
interply_first_layer = op.interply_first_layer;
interply_last_layer = op.interply_last_layer;
intraply_boundary = op.intraply_boundary;
%Define material indices of composite layers
layer1_i = fn_matl_i(matls, op.layer1);
layer2_i = fn_matl_i(matls, op.layer2);
interply_layer1_i = fn_matl_i(matls, op.interply_layer1);
interply_layer2_i = fn_matl_i(matls, op.interply_layer2);
intraply_layer1_i = fn_matl_i(matls, op.intraply_layer1);
intraply_layer2_i = fn_matl_i(matls, op.intraply_layer2);
%Validate input
if rem(n_ply_layers, 2)
    error("n_ply_layers must be even")
end
if op.upper_water_present
    error('fn_set_ply_material_v2 does not support op.upper_water_present, use v1')
end
if ~interply_boundary
    interply_first_layer = 0;
    interply_last_layer = 0;
end
if abs(interply_last_layer) > 1 || abs(interply_first_layer) > 1
    error('Interply_last/first_layer have to be either 0 or 1')
end
%Calculate ply layer heights so ply material is applied to specific layers
el_height = mod.el_height; %Calculated in fn_isometric_structured mesh
if op.interply_boundary
    n_interply_layers = (n_ply_layers - 1) + interply_last_layer + interply_first_layer;
else
    n_interply_layers = 0;
end
ply_height = (geom.specimen_height - n_interply_layers*el_height)/n_ply_layers;
height_offset = rem(ply_height, el_height);
%Alternate between upper and lower ply heights
ply_height_lower = ply_height - height_offset;
ply_height_upper = ply_height_lower + el_height;

% figure;
% %Display options
% display_options.interface_el_col = 'b';
% display_options.draw_elements = 0;
% %Plot geometry
% fn_show_geometry(mod, matls, display_options);
% hold on

%Assign materials
matl_oscillator = 1;
height_completed = 0;
%Loop over ply types
for ply_type = 1:n_ply_layers/n_plys_per_type
    %This changes ply material when needed
    if matl_oscillator == 1
        %first material
        matl_i = layer1_i;
        interply_matl_i = interply_layer1_i;
        intraply_matl_i = intraply_layer1_i;
    else
        %second material
        matl_i = layer2_i;
        interply_matl_i = interply_layer2_i;
        intraply_matl_i = intraply_layer2_i;
    end
    %Apply resin to first layer if enabled
    if interply_boundary && interply_first_layer && height_completed == 0
        mod = set_target_layer_material(mod, geom, interply_matl_i, el_height, height_completed);
        %Update height again
        height_completed = height_completed + el_height;
    end

    %Loop over layers in stack
    for layer_in_type = 1:n_plys_per_type
        target_layer = (ply_type - 1)*n_plys_per_type + layer_in_type;
        %Stay close to true ply layer height using upper and lower heights
        true_height = (target_layer - 1)*ply_height + (target_layer - 1 + interply_first_layer)*el_height;
        if height_completed > true_height
            ply_target_height = ply_height_lower;
        else
            ply_target_height = ply_height_upper;
        end

        % fn_show_geometry(mod, matls, display_options);
        % yline(model_height - true_height, 'Color','r')
        % yline(model_height - height_completed, 'Color','g')

        %Set materials in target layer
        mod = set_target_layer_material(mod, geom, matl_i, ply_target_height, height_completed);
        %Update height of material changed in specimen
        height_completed = height_completed + ply_target_height;
        
        %Set interply boundaries if enabled
        if interply_boundary
            if layer_in_type < n_plys_per_type || (ply_symmetry && (target_layer == n_ply_layers/2))
                if intraply_boundary
                    %Set boundary as intraply material
                    mod = set_target_layer_material(mod, geom, intraply_matl_i, el_height, height_completed);
                else
                    %Set boundary as intraply material
                    mod = set_target_layer_material(mod, geom, interply_matl_i, el_height, height_completed);
                end
                %Update height again
                height_completed = height_completed + el_height;
            elseif layer_in_type == n_plys_per_type && ~(interply_last_layer == 0 && target_layer == n_ply_layers) %weird logic but it works!
                %Set interlayer boundary
                mod = set_target_layer_material(mod, geom, interply_matl_i, el_height, height_completed);
                %Update height again
                height_completed = height_completed + el_height;
            end
        end
    end

    %Flip material type after n_plys_per_type loops completed
    matl_oscillator = -matl_oscillator;
    %CFRP specimen is symmetric
    if ply_symmetry && (target_layer == n_ply_layers/2)
        matl_oscillator = -matl_oscillator;
    end
end


% hold off

%% Calculate Model Accuracy Metrics

%Calculate deviation from true specimen height
model_specimen_height = height_completed;
model_accuracy = round(100*model_specimen_height/geom.specimen_height, 2);
fprintf("Specimen model height accuracy: %.2f%%\n", model_accuracy)
if abs(100 - model_accuracy) > 5
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('CAUTION: Specimen model accuracy is less than 95%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%Caluculate percentage of model made up of interply boundaries if enabled
if interply_boundary
    %Find number of elements used in specimen and number used for interply elements 
    layer_els = 0;
    for matl_i = [layer1_i, layer2_i]
        layer_els = layer_els + sum(mod.el_mat_i == matl_i);
    end
    interply_els = 0;
    for matl_i = [interply_layer1_i, interply_layer1_i]
        interply_els = interply_els + sum(mod.el_mat_i == matl_i);
    end
    total_els = layer_els + interply_els;

    %Calculate differences in rho and D by taking average change in rho and D
    rho_impact = mean([matls(interply_layer1_i).rho/matls(layer1_i).rho,...
        matls(interply_layer2_i).rho/matls(layer2_i).rho]);
    D_impact = mean([norm(matls(interply_layer1_i).D)/norm(matls(layer2_i).D),...
            norm(matls(interply_layer2_i).D)/norm(matls(layer2_i).D)]);
    %Print results
    interply_perc = round(100*interply_els/total_els,2);
    fprintf("Interply elements: Proportion:       %.2f%%\n", interply_perc)
    fprintf("                   Density impact:   %.2f%%\n", (rho_impact - 1)*interply_perc)
    fprintf("                   Stiffness impact: ~%.2f%%\n", (D_impact - 1)*interply_perc) %stiffness impact is approximate
end


%% Return top of specimen (because new height != defined specimen height)
top_of_specimen = height_completed;


end


function mod = set_target_layer_material(mod, geom, matl_i, target_height, height_completed)

%Extract key geometry consts
model_height = geom.model_height;
specimen_width = geom.specimen_width;
safety_margin = geom.safety_margin;

%Boundary points of target layer
bdry_pts = [0,               0
            specimen_width,  0
            specimen_width,  -target_height
            0,               -target_height];
%Translate upwards by height completed (starts from bottom of specimen)
bdry_pts = bdry_pts + (model_height - height_completed)*[0 1
                                                         0 1
                                                         0 1
                                                         0 1];
%Apply safety margin
bdry_pts = bdry_pts + safety_margin*[-1 0
                                      1 0
                                      1 0
                                     -1 0];
%Change materials
mod = fn_set_els_inside_bdry_to_mat(mod, bdry_pts, matl_i);
end