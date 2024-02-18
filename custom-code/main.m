clear;
close all;
restoredefaultpath;
addpath('../code');

%% Key modelling parameters

%Simulation Resolution
% --------------------------------------
%Elements per wavelength
els_per_wavelength = 5; %10 is default
time_step_safety_factor = 3; %3 is default
% --------------------------------------

%Location of transducer
src_displacement_from_surface = 0; %0 is default
%Specimen size and absorbing boundary + water regions as percent of specimen size
specimen_size = 4e-3; %[mm]
abs_bdry_thickness_perc = 0.05;
water_bdry_thickness_perc = 0.2;
%Material indices
ply_0_matl_i = 1;
ply_90_matl_i = 2;
water_matl_i = 3;
steel_matl_i = 4;
%Animate
animate = 0;

%% Define Materials

%Ply Material properties
E_fib = 163e+9; %[GPa]
G_fib = 5e+9;
v_fib = 0.32;
E_t = 12e+9;
G_t = 3.98e+9;
v_t = 0.024;

%Plys at 0 degrees orientation
ply_orientation = 0; %rotation of ply (0 or 90)
matls(ply_0_matl_i).rho = 1570; %Density
matls(ply_0_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation,E_fib,G_fib,v_fib,E_t,G_t,v_t);
matls(ply_0_matl_i).rayleigh_damping_coefs = [0 0]; %[alpha beta]
ply0col = hsv2rgb([2/3,0,0.80]);
matls(ply_0_matl_i).col = ply0col; %Colour for display
matls(ply_0_matl_i).name = 'Composite Ply Layer (0 deg)';
matls(ply_0_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Plys at 90 degrees orientation
ply_orientation = 90; %rotation of ply (0 or 90)
matls(ply_90_matl_i) = matls(ply_0_matl_i);
matls(ply_90_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation,E_fib,G_fib,v_fib,E_t,G_t,v_t);
matls(ply_90_matl_i).col = ply0col/2; %Colour for display
matls(ply_90_matl_i).name = 'Composite Ply Layer (90 deg)';

%Water
matls(water_matl_i).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = 1500 ^ 2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Water'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

%Steel (FOR DEBUGGING)
matls(steel_matl_i).rho = 8900; %Density
matls(steel_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_matl_i).name = 'Steel';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid


%% Define shape of model

%Define model parameters
water_brdy_thickness = water_bdry_thickness_perc*specimen_size;
model_size_w = specimen_size;
model_size_h = specimen_size + 2*water_brdy_thickness;
abs_bdry_thickness = abs_bdry_thickness_perc*specimen_size;

%Define size of model
bdry_pts = [
    0, 0 
    model_size_w, 0
    model_size_w, model_size_h
    0, model_size_h
    ];

%Define specimen size that will be water (water surrounds specimen)
wbt = water_brdy_thickness; %tmp for readability
specimen_brdy_pts = [
    0, wbt
    specimen_size, wbt
    specimen_size, wbt + specimen_size
    0, wbt + specimen_size];
%Define top of specimen for later use
top_of_specimen = specimen_brdy_pts(3,2);

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt, abt
    model_size_w - abt, abt
    model_size_w - abt, model_size_h - abt
    abt, model_size_h - abt];


%% Define Signal

src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)
centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;

%% Define elements and mesh

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% Define materials and absorbing layers

%First set material of all elements to water then set elements inside specimen 
%boundary to steel
mod.el_mat_i(:) = water_matl_i;
% mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, steel_matl_i);

% Set ply materials
[mod, top_of_specimen] = fn_set_ply_material(mod, ply_90_matl_i, ply_0_matl_i, specimen_brdy_pts);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

% Add porosity (WIP)
% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);


%% Define Probe (on mesh)

%Define a line along which sources will be placed to excite waves
src_end_pts = [
    0.25 * model_size_w, top_of_specimen - src_displacement_from_surface*wbt
    0.75 * model_size_w, top_of_specimen - src_displacement_from_surface*wbt];

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size, time_step_safety_factor);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%% DISPLAY MODEL

figure; 
display_options.interface_el_col = 'b';
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(mod, matls, display_options);

%% RUN THE MODEL

fe_options.field_output_every_n_frames = 10;
res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%% SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
if animate
    figure;
    display_options.interface_el_col = 'b';
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 1;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

