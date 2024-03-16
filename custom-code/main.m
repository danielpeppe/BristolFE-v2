clear;
close all;
% restoredefaultpath;
addpath('../code');

%% LOAD EXPERIMENTAL DATA

load('g4-s8_x-8000_z0_025.mat','exp_data');

%% DEFINE OUTPUT

op.porosity = 1; %[%]
%Dists
op.porosity_dist_sigma_tuner = 1.5;
op.porosity_dist_mu_tuner = 1.5; %0.5-1.5
%Morphology
op.porosity_r_min = 1e-6; %[m]
op.porosity_r_max = 3e-6;
%Matls
op.porosity_use_void = 0;
op.porosity_use_air = 0;
%   Density
op.porosity_use_density = 1;
op.porosity_sigma_tuner = 1;
op.porosity_mu_tuner = 1;
%   Damping
op.porosity_use_damping = 0;
%Output
% op_output.porosity_signal = 0;
op.porosity_plot_dists = 0;



%Resolution options
op.els_per_wavelength = 15; %increases are non-linear
op.time_step_safety_factor = 3; %12 if upper_water_present
anim_options.repeat_n_times = 10;
fe_options.field_output_every_n_frames = 50; %10 or inf is default (inf = no field output)
op.max_time = 3.5e-6; %3.5e-6;
%Output for each sim
op_output.justgeometry = 0; %disables other outputs
op_output.geometry = 1;
op_output.run_fea = 1;
op_output.plot_sim_data = 0;
op_output.plot_exp_data = 0;
op_output.animate = 0;
%Scale
op.plot_scale_dsps = 1.75; %1.75
op.plot_scale_time = 1;

%% TUNING

%Transducer options
op.aperture_n_els = 16; %number of elements
op.no_cycles = 3;
%Water options
op.water_rho_multiplier = 1.5;
op.water_D_multiplier = 1.5;
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
op.ply90_D_multiplier = 1.1;
op.ply0_D_multiplier = 1.1;
%   E_t Stiffness
op.ply90_E_t_multiplier = 1;
op.ply0_E_t_multiplier = 1;
%Damping options
op.rayleigh_quality_factor = 0.5; %inf disables damping (0.5 is light damping) (smaller = larger damping)
%Interply boundary options
op.interply_layer1 = 'resin';
op.interply_layer2 = 'resin';
op.interply_boundary = 1;
op.interply_el_thickness_perc = 0.05; %0.05
op.interply_rho_multiplier = 1;
op.interply_D_multiplier = 1;
op.interply_first_layer = 1;

%% PARAMS

%Define input parameters
op.params = [];
% op.params = [1 0];
% op.params = [10 12 15]; %els_per_wavelength
% op.params = [inf 0.1 0.5 1 5]; %damping
% op.params = {'resin','resin';'ply90','ply0'}; %intraply layers
% op.params = {[1,1],[1,2],[1,3],[2,1],[2,2],[2,3],[3,1],[3,2],[3,3]};
% op.params = {[0,90],[0,-90]};
% op.params = [0 1];
% op.params = [0.1 0.5 1 1.5 2.5];
% op.params = [1 0.1 0];
op.params = [3e-6 4e-6 5e-6];

%% RUNNING SIM

%Iterate sim for number of parameters
if isempty(op.params)
    fprintf("--------------------------- RUNNNING ONE SIM -----------------------------------\n")
    [op, op_output] = fn_set_options(op, op_output);
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
        p = op.params(i);
        %%%%%%%%%%%%%%%% Params start %%%%%%%%%%%%%%%%
        % op.els_per_wavelength = op.params(i);
        % op.intraply_layer1 = op.params{i,1}; op.intraply_layer2 = op.params{i,2};
        % op.ply0_orientation = op.params{1,i}(1); op.ply90_orientation = op.params{1,i}(2);
        op.porosity_r_max = p;
        %%%%%%%%%%%%%%%%% Params end %%%%%%%%%%%%%%%%%
        [op, op_output] = fn_set_options(op, op_output);
        fprintf("--------------------------------------------------------------------------------------\n")
        %Get results
        [res{1,i}, steps{1,i}] = run_sim(op, op_output, fe_options, anim_options, exp_data);
    end
    %Plot results
    fn_plot_signal(op, res, steps, exp_data)

    % if op_output.porosity_signal
    %     pristine_signal = res{1,1}{1}.dsps;
    %     porous_signal = res{1,2}{1}.dsps;
    %     porosity_signal{1,1}{1}.dsps = pristine_signal - porous_signal;
    %     fn_plot_signal(op, porosity_signal, steps, exp_data)
    % end
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
mat.ply90.rho = 1570 * op.ply90_rho_multiplier;
mat.ply90.D = op.ply90_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(op.ply90_orientation, E_fib, G_fib * op.ply90_G_x_multiplier, v_fib, E_t * op.ply90_E_t_multiplier, G_t, v_t);
mat.ply90.rayleigh_coefs = rayleigh_coefs;
mat.ply90.col = hsv2rgb([3/4,0.3,0.80]); %purple
mat.ply90.el_typ = 'CPE3';
%ply0
mat.ply0.rho = 1570 * op.ply0_rho_multiplier;
mat.ply0.D = op.ply0_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(op.ply0_orientation, E_fib, G_fib * op.ply0_G_x_multiplier, v_fib, E_t * op.ply0_E_t_multiplier, G_t, v_t);
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
mat.ply0b = mat.ply0;
mat.ply0b.rho = mat.ply0.rho * op.ply0b_interply_rho_multiplier;
mat.ply0b.D = mat.ply0.D * op.ply0b_interply_D_multiplier;
mat.ply0b.col = hsv2rgb([0,.75,.60]);
% ply90 boundary
mat.ply90b = mat.ply90;
mat.ply90b.rho = mat.ply90.rho * op.ply90b_interply_rho_multiplier;
mat.ply90b.D = mat.ply90.D * op.ply90b_interply_D_multiplier;
mat.ply90b.col = hsv2rgb([.40,.30,.60]);
%Plys not in between plys
mat.resin_intra = mat.resin;
mat.resin_intra.col = hsv2rgb([0.50,.75,.60]);
mat.resin_intra.rho = mat.resin.rho * op.intraply_rho_multiplier;
mat.resin_intra.D = mat.resin.D * op.intraply_D_multiplier;
%Water
% For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
% calcualted here from ultrasonic velocity (1500) and density (1000)
mat.water.rho = 1000 * op.water_rho_multiplier;
mat.water.D = 1500^2 * 1000 * op.water_D_multiplier;
mat.water.col = hsv2rgb([0.6,0.5,0.8]);
mat.water.el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid
%solid/fake water
if op.solidwater
    mat.solidwater.rho = 1000 * op.solidwater_rho_multiplier;
    mat.solidwater.D = op.solidwater_D_multiplier * fn_isotropic_plane_strain_stiffness_matrix(4e9, 0.3);
    mat.solidwater.col = hsv2rgb([0.6,0.75,0.8]);
    mat.solidwater.el_typ = 'CPE3';
end
%Steel
mat.steel.rho = 8900; %8900
mat.steel.D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
mat.steel.col = hsv2rgb([3/4,0.5,0.80]);
mat.steel.el_typ = 'CPE3';

%Define matls struct from mat
[matls, R_coefs] = fn_get_matls_struct(op, mat);

%% DEFINE SHAPE OF MODEL

%Define aperture size relative to exp probe
probe_width = op.scale_units * (exp_data.array.el_xc(end) - exp_data.array.el_xc(1));
aperture_width = double(probe_width * op.aperture_n_els/exp_data.num_els);
aperture_width_8 = double(probe_width * 8/exp_data.num_els);
%Define model parameters
water_brdy_thickness = op.water_bdry_thickness_perc * op.specimen_size;
if aperture_width > aperture_width_8
    model_width = (op.specimen_size + (aperture_width - aperture_width_8)) * op.model_width_multiplier;
else
    model_width = op.specimen_size * op.model_width_multiplier;
end
model_height = op.specimen_size + water_brdy_thickness * (op.upper_water_present + op.lower_water_present);
abs_bdry_thickness = op.abs_bdry_thickness_perc * op.specimen_size;

%Define size of model
model_bdry_pts = [
    0,            0 
    model_width,  0
    model_width,  model_height
    0,            model_height];

%Define specimen size that will be water (water surrounds specimen)
wbt = water_brdy_thickness; %tmp for readability
specimen_brdy_pts = [
    0,            wbt*op.lower_water_present
    model_width,  wbt*op.lower_water_present
    model_width,  wbt*op.lower_water_present + op.specimen_size
    0,            wbt*op.lower_water_present + op.specimen_size];

%Define top of specimen for later use
top_of_specimen = specimen_brdy_pts(3,2); 

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt,                abt*op.lower_water_present
    model_width - abt,  abt*op.lower_water_present
    model_width - abt,  model_height - op.upper_water_present*abt
    abt,                model_height - op.upper_water_present*abt];

%% DEFINE MESH

%Work out element size (slightly different from actual element size)
el_size = fn_get_suitable_el_size(matls, centre_freq, op.els_per_wavelength, op.scale_units);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(model_bdry_pts, el_size);

%% PRINT IMPORTANT INFORMATION

% fprintf('el_size: %.2f um\n', el_size * 1e6)
% [max_vel, min_vel] = fn_estimate_max_min_vels(matls, op.scale_units);
% avg_vel = (min_vel + max_vel) / 2;
% wavelengths = [min_vel max_vel avg_vel] / op.centre_freq * 1e6; %[um]
% fprintf('US wavelengths: min: %.2f, max: %.2f, avg: %.2f um\n', wavelengths)

%% SET MATERIALS IN MESH

%First set all elements to water
mod.el_mat_i(:) = fn_matl_i(matls,'water');
%Set upper elements to solid water if enabled
if op.solidwater
    solidwater_bdry = [0,            model_height/2
                       model_width,  model_height/2
                       model_width,  model_height
                       0,            model_height];
    mod = fn_set_els_inside_bdry_to_mat(mod, solidwater_bdry, fn_matl_i(matls,'solidwater'));
end
%Set specimen materials
if op.solid_specimen
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, fn_matl_i(matls,'ply90'));
elseif op.composite_specimen
    if op.upper_water_present
        [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, matls, specimen_brdy_pts);
    else
        %v2 does not suppot op.upper_water_present
        [mod, comp] = fn_set_ply_material_v2(mod, op, matls);
    end
end

%Add porosity
%Doesn't support intRAply layers or 2 types of intERply layers
[mod, matls] = fn_add_porosity_v4(op, mod, matls, comp);

%Add interface elements
mod = fn_add_fluid_solid_interface_els(mod, matls);
%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);


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
if strcmpi(op.src_matl,'solid_horizontal')
    src_end_pts = [
        0, 0.25 * model_height
        0, 0.75 * model_height];
else
    src_end_pts = [
        model_width/2 - aperture_width/2, top_of_specimen + src_offset
        model_width/2 + aperture_width/2, top_of_specimen + src_offset]; 
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
    % receiver_end_pts = src_end_pts - mod.el_height*[0 1; 0 1];
    receiver_end_pts = src_end_pts - op.specimen_size*[0 1; 0 1];
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
    % translate_time = 1.1964e-05; %Start of exp response

    % plot(steps{1}.load.time + translate_time * op.plot_scale_time, res_sum_dsps * op.plot_scale_dsps);
    plot(steps{1}.load.time, res_sum_dsps * op.plot_scale_dsps);
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
    scale_dsps = max(abs(res_sum_dsps))/max(abs(aperture_dsp_data)); %Scale exp response to match sim response
    translate_time = 1.1964e-05; %Start of exp response
    %Plot
    % plot(exp_data.time, aperture_dsp_data, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target');
    plot(exp_data.time - translate_time, aperture_dsp_data * scale_dsps, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target');

    hold off
    % figure; 
    % pulse = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);
    % plot(steps{1}.load.time + translate_time, pulse * 110)
    % hold on
    % plot(exp_data.time, aperture_dsp_data, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target')
end

%animate result
if animate
    figure;
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.norm_val = abs(median(res_sum_dsps)); %must be positive
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

end
