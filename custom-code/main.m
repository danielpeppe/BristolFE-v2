clear;
close all;
restoredefaultpath;
addpath('../code');

%% Key modelling parameters

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 6;

%% Define materials

%Material properties
steel_matl_i = 1;
matls(1).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(1).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
matls(1).rayleigh_damping_coefs = [0 0]; %[alpha beta]
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

water_matl_i = 2;
matls(water_matl_i).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = 1500 ^ 2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Water'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid


%% Define shape of model

%Define model parameters
specimen_size = 4e-3; %[mm]
water_brdy_thickness = 1e-3;
model_size_w = specimen_size;
model_size_h = specimen_size + 2*water_brdy_thickness;
abs_bdry_thickness = 1e-3;

%Define size of model
bdry_pts = [
    0, 0 
    model_size_w, 0
    model_size_w, model_size_h
    0, model_size_h
    ];

%Define specimen size that will be water (Water surrounds specimen)
wbt = water_brdy_thickness; %tmp for readability
specimen_brdy_pts = [
    0, wbt
    specimen_size, wbt
    specimen_size, wbt + specimen_size
    0, wbt + specimen_size];

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt, abt
    model_size_w - abt, abt
    model_size_w - abt, model_size_h - abt
    abt, model_size_h - abt];


%% Define Probe

%Define a line along which sources will be placed to excite waves
src_end_pts = [
    0.25 * model_size_w, model_size_h - wbt
    0.75 * model_size_w, model_size_h - wbt];

src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;

%% Define elements and mesh

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

% %DELETE ONCE BUG IS SOLVED (element size reduction for 9 to work)
% el_size = el_size * 0.85;

%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% Add porosity (WIP)

% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%% Defining materials and absorbing layers

%First set material of all elements to water then set elements inside specimen 
%boundary to steel
mod.el_mat_i(:) = water_matl_i;
mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, steel_matl_i);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);



%% Define Probe (on mesh)

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
figure; 
display_options.draw_elements = 1;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(mod, matls, display_options);

%--------------------------------------------------------------------------
%% RUN THE MODEL

fe_options.field_output_every_n_frames = 10;
res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%% SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
figure;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, display_options);
anim_options.repeat_n_times = 1;
fn_run_animation(h_patch, res{1}.fld, anim_options);

