clear;
close all;
restoredefaultpath;
addpath('../code');

%% LOAD EXPERIMENTAL DATA

load('g4-s8_x-8000_z0_025.mat');
%Signal
op.centre_freq = exp_data.array.centre_freq; %5e6 is default
no_cycles = 4; %4 is default

%% TUNING

%%%%%%%%%%%%% Tuning %%%%%%%%%%%%%
%Resolution
els_per_wavelength = 15; %10 is default (increases are non-linear)
time_step_safety_factor = 3; %3 is default (ensure reflections from 'unstable energy expansion' are avoided)
%Model options
op.model_size_w_multiplier = 1.5; %1.5 is default
op.abs_bdry_thickness_perc = 0.2; %0.2 is default (relative to specimen_size)
%Composite specimen options
op.layer1 = 'ply90';
op.layer2 = 'ply0';
%Ply options
op.n_plys_per_type = 2;
%   density
op.ply0_rho_multiplier = 1;
op.ply90_rho_multiplier = 1;
%   stiffness
op.ply90_D_multiplier = 1;
op.ply0_D_multiplier = 1;
op.ply90_shear_wave_velocity_by_E_t = 1; %1 is default
op.ply0_shear_wave_velocity_by_E_t = 1; %1 is default
%Damping options
op.rayleigh_quality_factor = inf; %inf disables damping
op.rayleigh_coefs = [0 1/(2*pi*op.centre_freq*(op.rayleigh_quality_factor * 1e4))]; %[alpha beta]
op.ply90_rayleigh_coefs = op.rayleigh_coefs;
op.ply0_rayleigh_coefs = op.rayleigh_coefs;
%Interply boundary options
op.interply_layer1 = 'resin';
op.interply_layer2 = 'resin';
op.interply_boundary = 1; %1 is default
op.interply_midway_boundary = 1; %1 is deafult
op.interply_every_layer = 1;
op.interply_rho_multiplier = 0.8;
op.interply_D_multiplier = 1;
%%%%%%%%%%%%% Tuning %%%%%%%%%%%%%

%% OTHER OPTIONS
%Paul demo
op.upper_water_present = 0;
op.water_interface_single = 0;
op.separate_transmitter = 0; %by 1 element
op.separate_receiver = 0;
%Sim options
fe_options.field_output_every_n_frames = 100; %10 or inf is default (inf = no field output)
max_time = 3.5e-6; %3.5e-6 is default (configures signal which in turn sets FEA time)
op.aperture_n_els = 8; %number of elements
%Output
geometry = 0; %disables other outputs
run_fea = 1;
plot_sim_data = 1;
plot_exp_data = 1;
animate = 1;

%Set default options and validate
op = fn_set_options(op);

%% DEFINE MATERIALS

%ply90
E_fib = 161e9; G_fib = 5.17e9; v_fib = 0.32; E_t = 11.38e9; G_t = 3.98e9; v_t = 0.436;
ply_orientation = 90; %rotation of ply (0 or 90 along z-axis)
mat.ply90.rho = 1570 * op.ply90_rho_multiplier;
mat.ply90.D = op.ply90_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t * op.ply90_shear_wave_velocity_by_E_t, G_t, v_t);
mat.ply90.rayleigh_coefs = op.ply90_rayleigh_coefs;
mat.ply90.col = hsv2rgb([3/4,0.3,0.80]); %purple
mat.ply90.el_typ = 'CPE3';
%ply0
ply_orientation = 0;
mat.ply0.rho = 1570 * op.ply0_rho_multiplier;
mat.ply0.D = op.ply0_D_multiplier * fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t * op.ply0_shear_wave_velocity_by_E_t, G_t, v_t);
mat.ply0.rayleigh_coefs = op.ply0_rayleigh_coefs;
mat.ply0.col = hsv2rgb([1/4,0,0.80]);
mat.ply0.el_typ = 'CPE3';
%Resin
mat.resin.rho = 1301 * op.interply_rho_multiplier;
mat.resin.D = op.interply_D_multiplier * fn_isotropic_plane_strain_stiffness_matrix(4.67e+9, 0.37); 
mat.resin.col = hsv2rgb([0,.75,.60]);
mat.resin.el_typ = 'CPE3';
%ply0 boundary
mat.ply0b = mat.ply0;
mat.ply0b.rho = mat.ply0.rho * op.interply_rho_multiplier;
mat.ply0b.D = mat.ply0.D * op.interply_D_multiplier;
mat.ply0b.col = hsv2rgb([0,.75,.60]);
%ply90 boundary
mat.ply90b = mat.ply90;
mat.ply90b.rho = mat.ply90.rho * op.interply_rho_multiplier;
mat.ply90b.D = mat.ply90.D * op.interply_D_multiplier;
mat.ply90b.col = hsv2rgb([.40,.30,.60]);
%Water
% For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
% calcualted here from ultrasonic velocity (1500) and density (1000)
mat.water.rho = 1000*op.water_rho_multiplier;
mat.water.D = 1500^2 * 1000;
mat.water.col = hsv2rgb([0.6,0.5,0.8]);
mat.water.el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid
%Steel (FOR DEBUGGING)
mat.steel.rho = 8900; %8900
mat.steel.D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
mat.steel.col = hsv2rgb([3/4,0.5,0.80]);
mat.steel.el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Get matls struct from mat struct
matls = fn_get_matls_struct(op,mat);

%% DEFINE SHAPE OF MODEL

%Define model parameters
water_brdy_thickness = op.water_bdry_thickness_perc*op.specimen_size;
model_size_w = op.specimen_size*op.model_size_w_multiplier;
model_size_h = op.specimen_size + water_brdy_thickness*(op.upper_water_present + op.lower_water_present);
abs_bdry_thickness = op.abs_bdry_thickness_perc*op.specimen_size;

%Define size of model
bdry_pts = [
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

%Work out element size
el_size = fn_get_suitable_el_size(matls, op.centre_freq, els_per_wavelength);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% DEFINE MATERIALS

%First set all elements to water
mod.el_mat_i(:) = fn_matl_i(matls,'water');
%Set specimen materials
if op.solid_specimen
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, fn_matl_i('steel'));
end
if op.composite_specimen
    [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, matls, specimen_brdy_pts);
end
%Add interface elements
mod = fn_add_fluid_solid_interface_els(mod, matls);
%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%% ADD POROSITY (WIP)
% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%% DEFINE PROBE END POINTS

%Define aperture
probe_width = exp_data.array.el_xc(end) - exp_data.array.el_xc(1);
aperture_width = probe_width * op.aperture_n_els/exp_data.num_els;
%Define a line along which sources will be placed to excite waves
if op.upper_water_present
    %Redefine src_dir for fluid
    src_dir = 4;
    %Redefine top of specimen
    top_of_specimen = new_top_of_specimen;
    %Adjust src offset if upper water is present
    if op.water_interface_single
        src_offset = mod.el_height;
    elseif op.water_interface_perc
        src_offset = op.water_interface_perc*wbt;
    else
        error('Critical Option Error!: problem with water options')
    end
else
    src_dir = 2;
    src_offset = 0;
end

%Define src end points
if ~op.horizontal_src
    src_end_pts = [
        model_size_w/2 - aperture_width/2, top_of_specimen + src_offset
        model_size_w/2 + aperture_width/2, top_of_specimen + src_offset];
else
    %Redefine src_dir for fluid
    src_dir = 1;
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
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size, time_step_safety_factor);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, op.centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
if op.separate_receiver
    receiver_end_pts = src_end_pts - mod.el_height*[0 1; 0 1];
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, receiver_end_pts(1, :), receiver_end_pts(2, :), el_size / 2);
else
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
end
steps{1}.mon.dfs = ones(size(steps{1}.mon.nds)) * src_dir;

%% DISPLAY MODEL

%Display options
display_options.interface_el_col = 'b';
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'b.';
%Plot geometry
if geometry
    figure; 
    h_patch = fn_show_geometry(mod, matls, display_options);
    return
end

%% RUN THE MODEL

if run_fea
    %Following relate to how absorbing regions are created by adding damping
    %matrix and reducing stiffness matrix to try and preserve acoustic impedance
    fe_options.damping_power_law = 3;
    fe_options.max_damping = 3.1415e+07;
    fe_options.max_stiffness_reduction = 0.01;
    
    %Run model
    res = fn_BristolFE_v2(mod, matls, steps, fe_options);
    res_sum_dsps = sum(res{1}.dsps); %tmp for readability

    %Measure instability using max, mean, median displacements (disabled)
    % fprintf('max dsp: %.2e, average dsp: %.2e median dsp: %.2e\n',max(res_sum_dsps), mean(res_sum_dsps), median(res_sum_dsps))
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
    aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx,aperture);
    aperture_time_data = sum(exp_data.time_data(:,aperture_data),2);
    %Scale dsp data
    scale_dsp = max(abs(res_sum_dsps))/max(abs(aperture_time_data)); %Scale exp response to match sim response
    translate_time = 1.194e-05; %Start of exp response
    %Plot
    plot(exp_data.time - translate_time, scale_dsp*aperture_time_data);
end

%animate result
if animate
    figure;
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 10;
    anim_options.norm_val = abs(median(res_sum_dsps)); %must be positive
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

