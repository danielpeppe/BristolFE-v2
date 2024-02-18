function mod = fn_set_ply_material(mod,ply_0_matl_i,ply_90_matl_i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%PROBLEMS:
%   1. size of element1a cannot be used because fn_isometric_structured_mesh
%      changes size of elements (and changes element numbers that are left
%      fairly arbitrarily.
%   2. mesh includes water, so number of elements in water must be
%      considered.
%   3. much easier to draw rectangle around elements.

% nodes_in_x_dir = mod.nds_in_x;
% nodes_in_y_dir = mod.nds_in_y;
el_mat_size = mod.el_mat_size;
els_in_x = el_mat_size(1);
els_in_y = el_mat_size(2);
els_in_mat = els_in_x*els_in_y;

%Calculate elements per ply layer as float first
els_per_ply_layer_float = els_in_y/32;
%Half as many element types per ply layer
els_type_per_ply_layer_float = els_per_ply_layer_float/2;

%Assign materials
y_checkpoint = 1;
matl_oscillator = 1;
for ply_types = 1:(32/4)

    if matl_oscillator > 0
        ply_matl_i = ply_0_matl_i;
    else
        ply_matl_i = ply_90_matl_i;
    end

    for ply_layer = 1:4
        if rem(ply_layer,2)
            els_type_per_ply_layer = floor(els_type_per_ply_layer_float);

            for x_checkpoint = 0:els_in_y:els_in_mat
                for element_target = y_checkpoint:(y_checkpoint + els_type_per_ply_layer)
                    %a type
                    mod.el_mat_i(x_checkpoint + element_target) = ply_matl_i;
                    %b type
                    mod.el_mat_i(els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                    %c type
                    mod.el_mat_i(2*els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                    %d type
                    mod.el_mat_i(3*els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                end
            end
        else
            els_type_per_ply_layer = ceil(els_type_per_ply_layer_float);

            for x_checkpoint = 0:els_in_y:els_in_mat
                for element_target = y_checkpoint:(y_checkpoint + els_type_per_ply_layer)
                    %a type
                    mod.el_mat_i(x_checkpoint + element_target) = ply_matl_i;
                    %b type
                    mod.el_mat_i(els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                    %c type
                    mod.el_mat_i(2*els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                    %d type
                    mod.el_mat_i(3*els_in_mat + x_checkpoint + element_target) = ply_matl_i;
                end
            end
        end 
        %Update starting point
        y_checkpoint = y_checkpoint + els_type_per_ply_layer;
    end
    matl_oscillator = matl_oscillator*-1;
end



end

