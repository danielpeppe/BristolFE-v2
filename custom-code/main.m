clear;
close all;
restoredefaultpath;
addpath('../code');

%% LOAD EXPERIMENTAL DATA

load('g4-s8_x-8000_z0_025.mat','exp_data');

%% DEFINE OUTPUT

%Resolution options
op.els_per_wavelength = 20; %increases are non-linear
op.time_step_safety_factor = 3; %12 if upper_water_present
anim_options.repeat_n_times = 10;
fe_options.field_output_every_n_frames = 50; %10 or inf is default (inf = no field output)
op.max_time = 3.5e-6;
%Output for each sim
op_output.justgeometry = 0; %disables other outputs
op_output.geometry = 1;
op_output.run_fea = 1;
op_output.plot_sim_data = 1;
op_output.plot_exp_data = 1;
op_output.animate = 1;

%% TUNING

%Model
op.model_size_w_multiplier = 1.5;
%Transducer options
op.aperture_n_els = 8; %number of elements
op.src_matl = 'solid';
%Water options
op.water_rho_multiplier = 1;
op.water_D_multiplier = 1;
%Composite specimen options
op.layer1 = 'ply90';
op.layer2 = 'ply0';
op.n_ply_layers = 32;
op.n_plys_per_type = 2;
op.ply_symmetry = 1;
%   Density
op.ply90_rho_multiplier = 1;
op.ply0_rho_multiplier = 1;
%   Stiffness
op.ply90_D_multiplier = 1;
op.ply0_D_multiplier = 1;
op.ply90_E_t_multiplier = 1;
op.ply0_E_t_multiplier = 1;
%Damping options
op.rayleigh_quality_factor = inf; %inf disables damping (0.5 is light damping) (smaller = larger damping)
%Interply boundary options
op.interply_layer1 = 'resin';
op.interply_layer2 = 'resin';
op.interply_boundary = 0;
op.interply_rho_multiplier = 1;
op.interply_D_multiplier = 1;
%Intraply bounday options
op.intraply_boundary = 0;
op.intraply_layer2 = 'resin_intra';
op.intraply_rho_multiplier = 1;
op.intraply_D_multiplier = 1;

%% RUNNING SIM

%Define input parameters
op.params = [];
% op.params = [1 0];
% op.params = [10 20 30]; %els_per_wavelength
% op.params = [inf 1]; %damping
% op.params = {'resin','resin';'ply90','ply0'}; %intraply layers
% op.params = [1 0.95 0.9 0.85 0.8]; %stiffness (D)
% op.params = {[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]};
% op.params = [0.9 0.95 1 1.05 1.1];
% op.params = [1 2]; %plys per type
% op.params = [1 4 8 12 14 16 32];

%Iterate sim for number of parameters
if isempty(op.params)
    fprintf("--------------------------- RUNNNING ONE SIM -----------------------------------\n")
    op = fn_set_options(op, op_output);
    fprintf("--------------------------------------------------------------------------------\n")
    
    %Get results
    [res{1}, steps{1}] = run_sim(op, op_output, fe_options, anim_options, exp_data);
else
    fprintf("--------------------------- RUNNNING MULTIPLE SIMS -----------------------------------\n")
    n = length(op.params);
    res = cell(1, n);
    steps = {1, n};
    for i = 1:n
        fprintf("----------------------------------- No %d/%d ------------------------------------------\n", i, n)
        %%%%%%%%%%%%%%%% Params start %%%%%%%%%%%%%%%%
        % op.els_per_wavelength = op.params(i);
        % op.rayleigh_quality_factor = op.params(i); %inf disables damping (0.5 is light damping)
        % op.interply_first_layer = op.params(i);
        % op.interply_boundary = op.params(i);
        % op.intraply_layer1 = op.params{i,1}; op.intraply_layer2 = op.params{i,2};
        % op.ply0_E_t_multiplier = op.params(i);
        % op.ply0_rho_multiplier = op.params(i);
        % op.ply90_G_x_multiplier = op.params(i);
        % op.ply0_G_x_multiplier = op.params(i);
        % op.n_plys_per_type = op.params(i);
        op.aperture_n_els = op.params(i);
        % op.ply_symmetry = op.params(i);
        %%%%%%%%%%%%%%%%% Params end %%%%%%%%%%%%%%%%%
        op = fn_set_options(op, op_output);
        fprintf("--------------------------------------------------------------------------------------\n")
        
        %Get results
        [res{1,i}, steps{1,i}] = run_sim(op, op_output, fe_options, anim_options, exp_data);
    end

    %Plot results
    fn_plot_signal(op, res, steps, exp_data, op.params)
end



function [res, steps] = run_sim(op, op_output, fe_options, anim_options, exp_data)
%% REDEFINE OPTIONS

%Paul Options
centre_freq = op.centre_freq;
no_cycles = op.no_cycles;
max_time = op.max_time;
%Options
geometry = op_output.geometry;
justgeometry = op_output.justgeometry;
run_fea = op_output.run_fea;
plot_sim_data = op_output.plot_sim_data;
plot_exp_data = op_output.plot_exp_data;
animate = op_output.animate;
if animate && fe_options.field_output_every_n_frames == inf || anim_options.repeat_n_times == 0
    error('fe_options.field_output_every_n_frames == inf, or, anim_options.repeat_n_times == 0, so animation will NOT be shown')
elseif ~animate
    fe_options.field_output_every_n_frames = inf;
end


%% DEFINE SPECIMEN MATERIALS

%Define rayleigh coefs
rayleigh_coefs = [0 1/(2*pi*centre_freq*(op.rayleigh_quality_factor * 1e4))]; %[alpha beta]

%ply90
E_fib = 161e9; G_fib = 5.17e9; v_fib = 0.32; E_t = 11.38e9; G_t = 3.98e9; v_t = 0.436;
ply_orientation = 90; %rotation of ply (0 or 90 along z-axis)
mat.ply90.rho = 1570 * op.ply90_rho_multiplier;
mat.ply90.D = op.ply90_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib * op.ply90_G_x_multiplier, v_fib, E_t * op.ply90_E_t_multiplier, G_t, v_t);
mat.ply90.rayleigh_coefs = rayleigh_coefs;
mat.ply90.col = hsv2rgb([3/4,0.3,0.80]); %purple
mat.ply90.el_typ = 'CPE3';
%ply0
ply_orientation = 0;
mat.ply0.rho = 1570 * op.ply0_rho_multiplier;
mat.ply0.D = op.ply0_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib * op.ply90_G_x_multiplier, v_fib, E_t * op.ply0_E_t_multiplier, G_t, v_t);
mat.ply0.rayleigh_coefs = rayleigh_coefs;
mat.ply0.col = hsv2rgb([1/4,0,0.80]);
mat.ply0.el_typ = 'CPE3';
%Resin
mat.resin.rho = 1301 * op.interply_rho_multiplier;
mat.resin.D = op.interply_D_multiplier * fn_isotropic_plane_strain_stiffness_matrix(4.67e+9, 0.37); 
mat.resin.col = hsv2rgb([0,.75,.60]);
mat.resin.el_typ = 'CPE3';

%% DEFINE BOUNDARY MATERIALS

%ply0 boundary
% mat.ply0b = mat.ply0;
% mat.ply0b.rho = mat.ply0.rho * op.interply_rho_multiplier;
% mat.ply0b.D = mat.ply0.D * op.interply_D_multiplier;
% mat.ply0b.col = hsv2rgb([0,.75,.60]);
%ply90 boundary
% mat.ply90b = mat.ply90;
% mat.ply90b.rho = mat.ply90.rho * op.interply_rho_multiplier;
% mat.ply90b.D = mat.ply90.D * op.interply_D_multiplier;
% mat.ply90b.col = hsv2rgb([.40,.30,.60]);
%Plys not in between plys
mat.resin_intra = mat.resin;
mat.resin_intra.col = hsv2rgb([0.50,.75,.60]);
mat.resin_intra.rho = mat.resin.rho * op.intraply_rho_multiplier;
mat.resin_intra.D = mat.resin.D * op.intraply_D_multiplier;
%Water
% For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
% calcualted here from ultrasonic velocity (1500) and density (1000)
mat.water.rho = 1000 * op.water_rho_multiplier;
mat.water.D = 1500^2 * 1000;
mat.water.col = hsv2rgb([0.6,0.5,0.8]);
mat.water.el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid
%solid/fake water
% if op.solidwater
%     mat.solidwater.rho = 1000 * op.solidwater_rho_multiplier;
%     solidwater_K = mat.water.D;
%     solidwater_v = 0.3;
%     mat.solidwater.D = op.solidwater_D_multiplier * fn_isotropic_plane_strain_stiffness_matrix(3*solidwater_K*(1-2*solidwater_v), solidwater_v);
%     mat.solidwater.col = hsv2rgb([0.6,0.75,0.8]);
%     mat.solidwater.el_typ = 'CPE3'; %AC2D3 must be the element type for a fluid
% end
%Steel
mat.steel.rho = 8900; %8900
mat.steel.D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
mat.steel.col = hsv2rgb([3/4,0.5,0.80]);
mat.steel.el_typ = 'CPE3';

%Get matls struct from mat struct
matls = fn_get_matls_struct(op, mat);

%% DEFINE SHAPE OF MODEL

%Define aperture size relative to exp probe
probe_width = op.scale_units * (exp_data.array.el_xc(end) - exp_data.array.el_xc(1));
aperture_width = double(probe_width * op.aperture_n_els/exp_data.num_els);
aperture_width_8 = double(probe_width * 8/exp_data.num_els);
%Define model parameters
water_brdy_thickness = op.water_bdry_thickness_perc * op.specimen_size;
if aperture_width > aperture_width_8
    model_size_w = (op.specimen_size + (aperture_width - aperture_width_8)) * op.model_size_w_multiplier;
else
    model_size_w = op.specimen_size * op.model_size_w_multiplier;
end
model_size_h = op.specimen_size + water_brdy_thickness * (op.upper_water_present + op.lower_water_present);
abs_bdry_thickness = op.abs_bdry_thickness_perc * op.specimen_size;

%Define size of model
model_bdry_pts = [
    0,            0 
    model_size_w, 0
    model_size_w, model_size_h
    0,            model_size_h];

%Define specimen size that will be water (water surrounds specimen)
wbt = water_brdy_thickness; %tmp for readability
specimen_brdy_pts = [
    0,            wbt*op.lower_water_present
    model_size_w, wbt*op.lower_water_present
    model_size_w, wbt*op.lower_water_present + op.specimen_size
    0,            wbt*op.lower_water_present + op.specimen_size];

%Define top of specimen for later use
top_of_specimen = specimen_brdy_pts(3,2);

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt,                abt*op.lower_water_present
    model_size_w - abt, abt*op.lower_water_present
    model_size_w - abt, model_size_h - op.upper_water_present*abt
    abt,                model_size_h - op.upper_water_present*abt];

%% DEFINE MESH

%Work out element size (slightly different from actual element size)
el_size = fn_get_suitable_el_size(matls, centre_freq, op.els_per_wavelength, op.scale_units);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(model_bdry_pts, el_size);

%% DEFINE MATERIALS

%First set all elements to water
mod.el_mat_i(:) = fn_matl_i(matls,'water');
%Set upper elements to solid water if enabled
if op.solidwater
    solidwater_bdry = [0,            model_size_h/2
                       model_size_w, model_size_h/2
                       model_size_w, model_size_h
                       0,            model_size_h];
    mod = fn_set_els_inside_bdry_to_mat(mod, solidwater_bdry, fn_matl_i(matls,'solidwater'));
end
%Set specimen materials
if op.solid_specimen
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, fn_matl_i(matls,'steel'));
elseif op.composite_specimen
    if op.upper_water_present
        [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, matls, specimen_brdy_pts);
    else
        %v2 does not suppot op.upper_water_present
        [mod, new_top_of_specimen] = fn_set_ply_material_v2(mod, op, matls, specimen_brdy_pts, model_size_h);
    end
end

%Add interface elements
mod = fn_add_fluid_solid_interface_els(mod, matls);
%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%% ADD POROSITY (WIP)

% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%% DEFINE PROBE END POINTS

%Define a line along which sources will be placed to excite waves
if op.upper_water_present
    %Redefine top of specimen
    top_of_specimen = new_top_of_specimen;
    %Adjust src offset if upper water is present
    if op.water_interface_single
        src_offset = mod.el_height;
    elseif op.water_interface_perc
        src_offset = op.water_interface_perc*wbt;
    else
        src_offset = 0;
    end
else
    src_offset = 0;
end

%Define src end points
if ~strcmpi(op.src_matl,'horizontal')
    src_end_pts = [
        model_size_w/2 - aperture_width/2, top_of_specimen + src_offset
        model_size_w/2 + aperture_width/2, top_of_specimen + src_offset];
else
    src_end_pts = [
        0, 0.25*model_size_h
        0, 0.75*model_size_h];
end

%% DEFINE LOADS

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
if op.separate_transmitter
    transmitter_end_pts = src_end_pts - mod.el_height*[0 1; 0 1];
    steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, transmitter_end_pts(1, :), transmitter_end_pts(2, :), el_size / 2);
else
    steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
end
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * op.src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size, op.scale_units, op.time_step_safety_factor);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
if op.separate_receiver
    receiver_end_pts = src_end_pts - mod.el_height*[0 1; 0 1];
    % receiver_end_pts = src_end_pts - 3.5e-3*[0 1; 0 1];
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, receiver_end_pts(1, :), receiver_end_pts(2, :), el_size / 2);
else
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
end
steps{1}.mon.dfs = ones(size(steps{1}.mon.nds)) * op.src_dir;

%% DISPLAY MODEL

%Display options
display_options.interface_el_col = 'b';
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'b.';
%Plot geometry
if justgeometry
    figure; 
    h_patch = fn_show_geometry(mod, matls, display_options);
    res{1} = {0};
    return
elseif geometry
    figure; 
    h_patch = fn_show_geometry(mod, matls, display_options);
end

%% RUN THE MODEL

if run_fea
    %Following relate to how absorbing regions are created by adding damping
    %matrix and reducing stiffness matrix to try and preserve acoustic impedance
    fe_options.damping_power_law = 3;
    fe_options.max_damping = 3.1415e+07 / op.scale_units;
    fe_options.max_stiffness_reduction = 0.01;
    %Run model
    res = fn_BristolFE_v2(mod, matls, steps, fe_options);
    res_sum_dsps = sum(res{1}.dsps); %tmp for readability
end

%% SHOW THE RESULTS

%Plot sim results
if plot_sim_data
    %Show the history output as a function of time - here we just sum over all 
    %the nodes where displacments were recorded
    figure;
    plot(steps{1}.load.time, res_sum_dsps);
    xlabel('Time (s)')
    ylabel('Magnitude (-)')
end
%Plot experimental data on top
if plot_exp_data
    hold on
    %Get time data for aperture
    aperture = 1:op.aperture_n_els;
    aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
    aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
    %Scale dsp data
    scale_dsp = max(abs(res_sum_dsps))/max(abs(aperture_dsp_data)); %Scale exp response to match sim response
    translate_time = 1.194e-05; %Start of exp response
    %Plot
    plot(exp_data.time - translate_time, scale_dsp*aperture_dsp_data, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target');
    hold off
end

%animate result
if animate
    figure;
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.norm_val = abs(median(res_sum_dsps)); %must be positive
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

end
