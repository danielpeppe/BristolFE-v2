clear;
close all;
restoredefaultpath;
addpath('../code');

%% LOAD EXPERIMENTAL DATA

load('g4-s8_x-8000_z0_025.mat');
%Signal
centre_freq = exp_data.array.centre_freq; %5e6 is default
no_cycles = 4; %4 is default

%% KEY MODELLING PARAMETERS

%%%%%%%%%%%%% Tuning %%%%%%%%%%%%%
%Resolution
els_per_wavelength = 5; %10 is default (increases are non-linear)
time_step_safety_factor = 3; %3 is default
%Model
op.model_size_w_multiplier = 1.5; %1.5 is default
op.abs_bdry_thickness_perc = 0.2; %0.2 is default (relative to specimen_size)
%Ply options
op.shear_wave_velocity_by_E_t = 1.17; %1 is default (adjusts E_t stiffness) (1.27)
op.back_wall_reflection_by_water_density = 1; %1 is default
op.rayleigh_quality_factor = 3000;
op.rayleigh_coefs = [0 1/(2*pi*centre_freq*op.rayleigh_quality_factor)]; %[alpha beta]
%op.rayleigh_coefs = [0 0];
%Interply boundary
op.interply_boundary = 0;
op.interply_density_multiplier = 1.05;
op.interply_stiffness_multiplier = 1.05;
%%%%%%%%%%%%% Tuning %%%%%%%%%%%%%

%Active
% op.upper_water_present = 1;
% op.water_interface_single = 1;
op.n_plys_per_type = 2;

%Signal
op.separate_transmitter = 0;
op.separate_receiver = 0;
op.aperture_n_els = 8; %number of elements
%Material indices (ply Material layers = 1 and 2) (solid mateial = 3) (indices cannot be skipped)
ply90_matl_i = 1;
ply90boundary_matl_i = 2; %must be layer index + 1
ply0_matl_i = 3;
ply0boundary_matl_i = 4;
steel_matl_i = 5; %DEBUGGING
water_matl_i = 6;
%FEA options
fe_options.field_output_every_n_frames = 100; %10 or inf is default (inf = no field output)
max_time = 3.5e-6; %5e-6 is default (configures signal which in turn sets FEA time)
%Output
op.geometry = 0; %disables other outputs
op.run_fea = 1;
op.plot_sim_data = 1;
op.plot_exp_data = 1;
op.animate = 1;

%Set default options and validate
op = fn_set_options(op);

%% DEFINE MATERIAL

%Ply Material properties
% E_fib = 164e9; G_fib = 5e9; v_fib = 0.32; G_t = 3.98e9; v_t = 0.024;
% E_t = 12e9 * op.shear_wave_velocity_by_E_t; %wave velocity can be adjusted by transverse stiffness
E_fib = 86.65e9; G_fib = 5.17e9; v_fib = 0.435; G_t = 4.50e9; v_t = 0.042;
E_t = 13.44e9 * op.shear_wave_velocity_by_E_t; %wave velocity can be adjusted by transverse stiffness

%Plys at 0 degrees orientation (along z-axis)
ply_orientation = 0;
matls(ply0_matl_i).rho = 1570;
matls(ply0_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t, G_t, v_t);
matls(ply0_matl_i).rayleigh_coefs = op.rayleigh_coefs;
ply0col = hsv2rgb([3/4,0,0.80]);
matls(ply0_matl_i).col = ply0col;
matls(ply0_matl_i).name = 'Ply Layer (0 deg)';
matls(ply0_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Plys at 90 degrees orientation (along x-axis)
ply_orientation = 90; %rotation of ply (0 or 90)
matls(ply90_matl_i) = matls(ply0_matl_i);
matls(ply90_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t, G_t, v_t);
matls(ply90_matl_i).col = ply0col/2;
matls(ply90_matl_i).name = 'Ply Layer (90 deg)';

%Ply boundary materials
matls(ply0boundary_matl_i) = matls(ply0_matl_i);
matls(ply0boundary_matl_i).rho = matls(ply0_matl_i).rho * op.interply_density_multiplier;
matls(ply0boundary_matl_i).D = matls(ply0_matl_i).D * op.interply_stiffness_multiplier;
matls(ply0boundary_matl_i).col = hsv2rgb([0,.75,.60]);
matls(ply90boundary_matl_i) = matls(ply90_matl_i);
matls(ply90boundary_matl_i).rho = matls(ply90_matl_i).rho * op.interply_density_multiplier;
matls(ply90boundary_matl_i).D = matls(ply90_matl_i).D * op.interply_stiffness_multiplier;
matls(ply90boundary_matl_i).col = hsv2rgb([.40,.30,.60]);

%Water
matls(water_matl_i).rho = 1000*op.back_wall_reflection_by_water_density;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = 1500^2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Water'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

%Steel (FOR DEBUGGING)
matls(steel_matl_i).rho = 8900; %8900
matls(steel_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
matls(steel_matl_i).col = hsv2rgb([3/4,0.5,0.80]);
matls(steel_matl_i).name = 'Steel';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Resin (FOR DEBUGGING)
% matls(5).rho = 1301; %Density
% matls(5).D = fn_isotropic_plane_strain_stiffness_matrix(4.67e+9, v_t); 
% matls(5).col = hsv2rgb([2/3,0,0.80]); %Colour for display
% matls(5).name = 'Resin';
% matls(5).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

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
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% DEFINE MATERIALS

%First set all elements to water
mod.el_mat_i(:) = water_matl_i;
%Set specimen materials
if op.solid_specimen
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, 5);
end
if op.composite_specimen
    %Set ply materials (input 1 = layer 1, input 2 = layer 2)
    [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, 1, 3, specimen_brdy_pts);
end
%Add interface elements
mod = fn_add_fluid_solid_interface_els(mod, matls);
%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%% ADD POROSITY (WIP)
% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%% DEFINE PROBE END POINTS

%Define exp probe aperture
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
        error('Critical Option Error: problem if water options')
    end
else
    src_dir = 2;
    src_offset = 0;
    %src_offset = -0.5*model_size_h + wbt;
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
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

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

if op.geometry
    figure; 
    h_patch = fn_show_geometry(mod, matls, display_options);
end

%% RUN THE MODEL

if op.run_fea
    %Following relate to how absorbing regions are created by adding damping
    %matrix and reducing stiffness matrix to try and preserve acoustic impedance
    fe_options.damping_power_law = 3;
    fe_options.max_damping = 3.1415e+07;
    fe_options.max_stiffness_reduction = 0.01;
    
    %Run model
    res = fn_BristolFE_v2(mod, matls, steps, fe_options);
    res_sum_dsps = sum(res{1}.dsps); %save displacements in variable for readability

    %Measure instability using max, mean, median displacements (disabled)
    % fprintf('max dsp: %.2e, average dsp: %.2e median dsp: %.2e\n',max(res_sum_dsps), mean(res_sum_dsps), median(res_sum_dsps))
end

%% SHOW THE RESULTS

%Plot sim results
if op.plot_sim_data
    %Show the history output as a function of time - here we just sum over all 
    %the nodes where displacments were recorded
    figure;
    plot(steps{1}.load.time, res_sum_dsps);
    xlabel('Time (s)')
    ylabel('Magnitude (-)')
end
%Plot experimental data on top
if op.plot_exp_data
    hold on
    %Get time data for aperture
    aperture = 1:op.aperture_n_els;
    aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx,aperture);
    aperture_time_data = sum(exp_data.time_data(:,aperture_data),2);
    %Scale dsp data
    scale_dsp = max(res_sum_dsps)/max(aperture_time_data); %Scale exp response to match sim response
    translate_time = 1.194e-05; %Start of exp response
    %Plot
    plot(exp_data.time - translate_time, scale_dsp*aperture_time_data);
end

%animate result
if op.animate
    figure;
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 10;
    anim_options.norm_val = abs(median(res_sum_dsps)); %must be positive
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

