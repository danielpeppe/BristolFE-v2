function [mod, matls, op] = fn_add_porosity_v4(mod, op, matls, comp)
%FN_GEN_VOID Summary of this function goes here
%   Detailed explanation goes here

%Validate input
porosity = op.porosity;
porosity_r_min = op.porosity_r_min;
porosity_r_max = op.porosity_r_max;
if porosity == 0
    return
elseif porosity_r_min > porosity_r_max
    error('op.porosity_r_min > op.porosity_r_max is required')
else
    porosity = porosity/100; %perc -> decimal
end
fprintf("Add Porosity (v4)")


%Get material indices from comp struct (defined in fn_set_ply_material)
layer1_i = comp.matl_i.layer1_i;
layer2_i = comp.matl_i.layer2_i;
interply_layer_i = comp.matl_i.interply_layer1_i;
ply_width = comp.specimen_width;
matls_for_pores_arr = fn_matl_i(matls, op.porosity_matls_for_pores);

%% DEFINE PORE MATERIALS

if op.porosity_use_void
    fprintf(" using voids\n")
elseif op.porosity_use_density

    %Print important info
    porosity_r_absolute_max = sqrt(mod.el_height^2/(sqrt(3)*pi));
    if porosity_r_max > porosity_r_absolute_max
        error('porosity_r_max is too high. Decrease porosity_r_max.')
    end
    % porosity_r_range = [porosity_r_min, porosity_r_max, porosity_r_absolute_max] * 1e6;
    % fprintf("\n")
    % fprintf("using density: rad range: [%.2f, %.2f] (max: %.2f)\n", porosity_r_range)
    % fprintf("               el_height: %.2f\n", mod.el_height * 1e6)
    % fprintf("               max rho reduction: %.2f times\n", 1/(1 - pi*sqrt(3)*((porosity_r_max^2/mod.el_height^2))))
    

    %Loop over pore types
    col_sat_arr = linspace(0.25, 1, op.porosity_n_pore_matls);
    %col_brightness_arr = linspace(0.5, 0.7, numel(op.porosity_matls_for_pores));
    pore_r_arr = linspace(porosity_r_min, porosity_r_max, op.porosity_n_pore_matls);
    mat = struct();
    for i = 1:op.porosity_n_pore_matls
        %Loop over materials which will contrain porosity
        for ii = 1:length(op.porosity_matls_for_pores)
            ply_mat_i = fn_matl_i(matls, op.porosity_matls_for_pores(ii));
            pore_mat_field = strcat(op.porosity_matls_for_pores(ii),"_pore",string(i));
            %Define scale factor in terms of radius of pore
            scale_factor = 1 - pi*sqrt(3)*((pore_r_arr(i)^2/mod.el_height^2));
            %Define density
            mat.(pore_mat_field).rho = matls(ply_mat_i).rho * scale_factor;
            % mat.(pore_mat_field).rho = matls(ply_mat_i).rho;
            %Define stiffness as proportional to drop in density (to preserve impedance)
            mat.(pore_mat_field).D = matls(ply_mat_i).D * scale_factor;
            % mat.(pore_mat_field).D = matls(ply_mat_i).D;
            %Define everything else
            mat.(pore_mat_field).rayleigh_coefs = matls(ply_mat_i).rayleigh_coefs / scale_factor * op.porosity_damping_tuner ;
            mat.(pore_mat_field).col = hsv2rgb([.169, col_sat_arr(i), 1]); %just make sure colours are distinct, strange calculation here
            mat.(pore_mat_field).el_typ = 'CPE3';
        end
    end
    %Define matls struct from mat
    [matls, Z, R_coefs] = fn_get_matls_struct(op, mat, matls);
end

%% NUMBER OF PORES

%Pore radius is linearly proportional to level of porosity within range
porosity_range = op.porosity_range(2) - op.porosity_range(1);
mu1 = porosity_r_min + (porosity_r_max - porosity_r_min) * (op.porosity/porosity_range);
sigma1 = (porosity_r_max - porosity_r_min)/4 * op.porosity_r_sigma_tuner;
pore_r_pd = makedist('Normal', 'mu', mu1, 'sigma', sigma1);
% Truncate the normal distribution
pore_r_trunc_pd = truncate(pore_r_pd, porosity_r_min, porosity_r_max);

%Define composite density as weight ratio of ply0 and resin density
ply_rho = matls(layer1_i).rho; %can be either ply0 or ply90, doesn't matter which
resin_rho = matls(interply_layer_i).rho;
comp.pristine_rho = comp.interply_volume_frac*resin_rho + comp.ply_volume_frac*ply_rho;
%Turn porosity value into an equivalent decrease in ply0 layer density
true_porous_ply_rho = (comp.pristine_rho*(1 - porosity) - resin_rho*comp.interply_volume_frac)/comp.ply_volume_frac;

%Calculate number of ply elements
n_ply_els = 0;
for matl_i = [layer1_i layer2_i]
    n_ply_els = n_ply_els + sum(mod.el_mat_i == matl_i);
end

%Calculate numer of pore elements
if op.porosity_use_void
    total_n_pores = round(n_ply_els * (ply_rho - true_porous_ply_rho)/ply_rho);
    n_ply_els = n_ply_els - total_n_pores;
elseif op.porosity_use_density
    %Generate random pore radii
    pore_r_rand = random(pore_r_trunc_pd, n_ply_els, 1); %number of pores cannot exceed n_ply_els
    %Convert into densities
    el_rho_rand = ply_rho*(1 - pi*sqrt(3)*(((pore_r_rand/2).^2/mod.el_height^2))); %THIS DOESNT CONSIDER DIFFERENT PLY DENSITIES
    total_n_pores = 0;
    porous_ply_rho = ply_rho;
    %Generate list of random pore element densities to assign later and count number of pores required
    while porous_ply_rho > true_porous_ply_rho
        total_n_pores = total_n_pores + 1;
        pore_rho = el_rho_rand(total_n_pores);
        porous_ply_rho = (ply_rho*(n_ply_els - total_n_pores) + total_n_pores*pore_rho)/n_ply_els;
        n_ply_els = n_ply_els - 1;
    end
    pore_rho_arr = el_rho_rand(1:total_n_pores);

    %Print true porosity
    comp.porosity_rho = comp.interply_volume_frac*resin_rho + comp.ply_volume_frac*porous_ply_rho;
    actual_porosity = 100*(1 - comp.porosity_rho/comp.pristine_rho);
    fprintf('               Actual Porosity: %.2f%% w/ %d Pores (target: %.2f%%)\n', actual_porosity, total_n_pores, 100*porosity)
end

%% PORE SPATIAL DISTRIBUTION

%Create ply layer probability distribution
upper_layer = comp.ply_location_tracker{1, 1};
lower_layer = comp.ply_location_tracker{op.n_ply_layers, 2};
mu2 = (upper_layer + lower_layer)/2 * op.porosity_dist_mu_tuner;
sigma2 = (upper_layer - lower_layer)/4 * op.porosity_dist_sigma_tuner;
pore_dist_pd = makedist('Normal', 'mu', mu2, 'sigma', sigma2);
%Truncate the normal distribution
pore_dist_trunc_pd = truncate(pore_dist_pd, lower_layer, upper_layer);

%Assign pore material
% While loop is required since it is not guaranteed that the random pore
% assignment will result in all the required pores being placed, since only
% n_rand_locs pores are placed, many of which may be placed in invalid
% locations (in resin)
n_samples = round(total_n_pores * op.porosity_dist_n_samples_sf);
pore_el_i_valid = zeros(total_n_pores, 1);
while ismember(pore_el_i_valid, 0)

    %Get location where 'porous element' will be locates
    rand_height_arr = random(pore_dist_trunc_pd, n_samples, 1);
    rand_width_arr = ply_width * rand(n_samples, 1);
    %Build KD-tree with element centers
    kdTree = createns(mod.el_centres, 'NSMethod','kdtree');
    %Get random location within ply to place porosity
    rand_loc_arr = [rand_width_arr, rand_height_arr];
    %Bulk nearest neighbor search
    pore_el_i_arr = knnsearch(kdTree, rand_loc_arr, 'K', 1);

    %Remove duplicates
    pore_el_i_arr_un = unique(pore_el_i_arr);
    %Get chosen pore element materials
    pore_el_matl_i = mod.el_mat_i(pore_el_i_arr_un);
    pore_el_matl_i_is_pore_matl = ismember(pore_el_matl_i, matls_for_pores_arr);
    if sum(pore_el_matl_i_is_pore_matl) > total_n_pores
        %Set ply elements to pore materials
        pore_el_i_valid = pore_el_i_arr_un(pore_el_matl_i_is_pore_matl);
        shuffle_indices = randperm(numel(pore_el_i_valid), total_n_pores);
        pore_el_i_valid = pore_el_i_valid(shuffle_indices);
    else
        %Reset valid element indices
        pore_el_i_valid = zeros(total_n_pores, 1);
        fprintf("trying again...")
    end
end

%% ASSIGN PORE MATERIAL

if op.porosity_use_void
    %Remove pores from mesh
    mod.els(pore_el_i_valid, :) = [];
    mod.el_mat_i(pore_el_i_valid) = [];
    mod.el_centres(pore_el_i_valid, :) = [];
    [mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);
else
    %Extract pore materials from matls struct
    matls_names = {matls.name};
    for ply_mat = op.porosity_matls_for_pores
        ply_mat_c = char(ply_mat);
        ply_mat_i = [ply_mat_c,'_i'];
        pore_matls.(ply_mat_i) = 1;
        for i = 1:numel(matls)
            pore_mat_num = pore_matls.(ply_mat_i);
            num_cutoff = floor(abs(log10(pore_mat_num))) + 1;
            if strcmp([ply_mat_c,'_pore'], matls_names{i}(1:end - num_cutoff))
                pore_matls.(ply_mat)(pore_mat_num) = matls(i);
                pore_matls.(ply_mat_i) = pore_mat_num + 1;
            end
        end
    end

    %Assign pore materials
    if op.porosity_use_density
        for i = 1:length(pore_el_i_valid)
            %Extract pore matl and distribution
            pore_rho = pore_rho_arr(i);
            pore_el_i = pore_el_i_valid(i);
            %Get correct ply material
            ply_mat = matls(mod.el_mat_i(pore_el_i)).name;
            %Find best suited pore material to assign to pore element
            [~, pore_mat_i_in_pore_matls] = min(abs([pore_matls.(ply_mat)(:).rho] - pore_rho));
            %Assign pore material
            pore_mat_i_in_matls = find(strcmp({matls.name}, pore_matls.(ply_mat)(pore_mat_i_in_pore_matls).name));
            mod.el_mat_i(pore_el_i) = pore_mat_i_in_matls;
        end
    end
    %Validate porosity material has been assigned to model elements
    if ismember(mod.el_mat_i, pore_el_i_arr_un)
        error('Porosity has not been added to model for some reason')
    end
end

%% RETURN POROSITY DATA

op.actual_porosity = actual_porosity;
op.total_n_pores = total_n_pores;

%% PLOT DISTRIBUTIONS

if op.porosity_plot_dists
    figure;
    % First subplot
    subplot(2,2,1);
    histogram(pore_r_rand, 'Normalization', 'pdf', 'BinWidth', mean(abs(diff(pore_r_rand)))/4)
    title('Pore Radius Distribution');
    % Second subplot
    subplot(2,2,2);
    histogram(el_rho_rand, 'Normalization', 'pdf', 'BinWidth', mean(abs(diff(el_rho_rand)))/4)
    title('Element Density Distribution');
    % Third subplot
    subplot(2,2,3);
    histogram(pore_rho_arr, 'Normalization', 'pdf', 'BinWidth', mean(abs(diff(pore_rho_arr)))/4)
    title('Actual Element Density Distribution');
    % Fourth subplot
    subplot(2,2,4);
    pore_dist_trunc_pd_rand = random(pore_dist_trunc_pd, 1, n_ply_els);
    histogram(pore_dist_trunc_pd_rand, 'Normalization', 'pdf', 'BinWidth', mean(abs(diff(pore_dist_trunc_pd_rand)))/4)
    title('Actual Pore Spatial Distribution');
end

fprintf("...completed\n")

end


