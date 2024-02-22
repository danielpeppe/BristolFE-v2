clear;
close all;
restoredefaultpath;
addpath('../code');

%% KEY MODELLING PARAMETERS

%Simulation Resolution --------------------------------------
els_per_wavelength = 10; %10 is default (increases are non-linear)
time_step_safety_factor = 12; %3 is default
%Simulation Parameters --------------------------------------
%Specimen size and absorbing boundary + water regions as percent of specimen size
specimen_size = 4e-3; %[mm]
abs_bdry_thickness_perc = 0.1; %0.1 is default (relative to water_boundary_size)
water_bdry_thickness_perc = 0.2; %0.2 is default (>abs_bdry_thickness_perc)
%Material indices (Ply Material layers = 1 and 2) (indices cannot be skipped)
ply90_matl_i = 1;
ply0_matl_i = 2;
water_matl_i = 3;
steel_matl_i = 4; %DEBUGGING
%FEA options
field_output_every_n_frames = 100; %10 or inf is default (inf = no field output)
use_gpu_if_present = 1;
%Switches ---------------------------------------------------
%Location of transducer
sw.water_interface_single = 0; %0 is default (1 separates transducer from specimen by 1 element)
sw.water_interface = 0.1; %0 is default (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
%Input
sw.steel_and_water = 0;
sw.plys_and_water = 1;
%Output
sw.geometry = 0;
sw.plot_exp_data = 1;
sw.run_fea = 1;
sw.animate = 1;
% -----------------------------------------------------------

%% CALCULATE MORE PARAMETERS

if sw.water_interface && sw.water_interface_single
    error('Choose either sw.water_interface OR sw.water_interface_single (not both)')
elseif sw.water_interface || sw.water_interface_single
    sw.water_present = 1;
else
    sw.water_present = 0;
end

%% DEFINE MATERIAL

%Steel (FOR DEBUGGING)
matls(steel_matl_i).rho = 8900;
matls(steel_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
matls(steel_matl_i).col = hsv2rgb([3/4,0.5,0.80]);
matls(steel_matl_i).name = 'Steel';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Ply Material properties
E_fib = 163e+09; %[GPa]
G_fib = 5e+09;
v_fib = 0.32;
E_t = 12e+09;
G_t = 3.98e+09;
v_t = 0.024;
%v_t = 0.4;
ply_rho = 1570; %[kg/m3]

%Plys at 0 degrees orientation (along z-axis)
ply_orientation = 0; %rotation of ply (0 or 90)
matls(ply0_matl_i).rho = ply_rho;
matls(ply0_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t, G_t, v_t);
%matls(ply0_matl_i).D = matls(steel_matl_i).D;
% matls(ply0_matl_i).rayleigh_damping_coefs = [0 0]; %[alpha beta]
ply0col = hsv2rgb([3/4,0,0.80]);
matls(ply0_matl_i).col = ply0col;
matls(ply0_matl_i).name = 'Ply Layer (0 deg)';
matls(ply0_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Plys at 90 degrees orientation (along x-axis)
ply_orientation = 90; %rotation of ply (0 or 90)
matls(ply90_matl_i) = matls(ply0_matl_i);
matls(ply90_matl_i).D = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation, E_fib, G_fib, v_fib, E_t, G_t, v_t);
%matls(ply90_matl_i).D = matls(steel_matl_i).D; 
matls(ply90_matl_i).col = ply0col/2;
matls(ply90_matl_i).name = 'Ply Layer (90 deg)';

%Water
matls(water_matl_i).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = 1500^2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Water'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid


%Resin (FOR DEBUGGING)
% matls(5).rho = 1301; %Density
% matls(5).D = fn_isotropic_plane_strain_stiffness_matrix(4.67e+9, v_t); 
% matls(5).col = hsv2rgb([2/3,0,0.80]); %Colour for display
% matls(5).name = 'Resin';
% matls(5).el_typ = 'CPE3'; %CPE3 must be the element type for a solid


%% DEFINE SHAPE OF MODEL

%Define model parameters
water_brdy_thickness = water_bdry_thickness_perc*specimen_size;
model_size_w = specimen_size;
if sw.water_present
    model_size_h = specimen_size + 2*water_brdy_thickness;
else
    model_size_h = specimen_size + water_brdy_thickness;
end
abs_bdry_thickness = abs_bdry_thickness_perc*specimen_size;

%Define size of model
bdry_pts = [
    0,            0 
    model_size_w, 0
    model_size_w, model_size_h
    0,            model_size_h];

%Define specimen size that will be water (water surrounds specimen)
wbt = water_brdy_thickness; %tmp for readability
specimen_brdy_pts = [
    0,            wbt
    model_size_w, wbt
    model_size_w, wbt + specimen_size
    0,            wbt + specimen_size];
%Define top of specimen for later use
top_of_specimen = specimen_brdy_pts(3,2);

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt,                abt
    model_size_w - abt, abt
    model_size_w - abt, model_size_h - sw.water_present*abt
    abt,                model_size_h - sw.water_present*abt];

%% DEFINE SIGNAL

centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;
if sw.water_interface || sw.water_interface_single
    src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)
else
    src_dir = 2;
end

%% DEFINE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% DEFINE MATERIALS AND POROSITY

mod.el_mat_i(:) = water_matl_i;
if sw.steel_and_water
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, steel_matl_i);
end
if sw.plys_and_water
    %Set ply materials (input 1 = layer 1, input 2 = layer 2)
    [mod, new_top_of_specimen] = fn_set_ply_material(mod, 1, 2, specimen_brdy_pts, 32, 4);
    %Redefine top of specimen
    top_of_specimen = new_top_of_specimen;
end
%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Add porosity (WIP)
% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);


%% DEFINE TRANSDUCER

%Define a line along which sources will be placed to excite waves
if sw.water_interface_single
    src_offset = mod.el_height;
elseif sw.water_interface
    src_offset = sw.water_interface*wbt;
else
    src_offset = 0;
    %src_offset = -0.5*model_size_h + wbt;
end

src_end_pts = [
    0.25 * model_size_w, top_of_specimen + src_offset
    0.75 * model_size_w, top_of_specimen + src_offset];

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

if sw.geometry
    figure; 
    display_options.interface_el_col = 'b';
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
end

%% RUN THE MODEL

if sw.run_fea
    %Define FE options
    fe_options.field_output_every_n_frames = field_output_every_n_frames;
    fe_options.use_gpu_if_present = use_gpu_if_present;
    
    fe_options.dof_to_use = []; %Blank uses all of available ones for all elements, subset can be set e.g. as [1,2]
    %Following relate to how absorbing regions are created by adding damping
    %matrix and reducing stiffness matrix to try and preserve acoustic impedance
    fe_options.damping_power_law = 3;
    fe_options.max_damping = 3.1415e+07;
    fe_options.max_stiffness_reduction = 0.01;
    %Run model
    res = fn_BristolFE_v2(mod, matls, steps, fe_options);

    fprintf('max dsp: %2f, average dsp: %2f\n',max(sum(res{1}.dsps)), mean(sum(res{1}.dsps)))
end

%% SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')
hold on

if sw.plot_exp_data
    load('g4-s8_x-8000_z0_025.mat');
    aperture = 1:8;
    i = ismember(exp_data.tx, aperture) & ismember(exp_data.rx,aperture);
    time_data_i = sum(exp_data.time_data(:,i),2);
    scale_dsp = max(sum(res{1}.dsps))/max(time_data_i);
    translate_time = 1.194e-05;
    plot(exp_data.time - translate_time, scale_dsp*time_data_i);
end
hold off

%sw.animate result
if sw.animate
    figure;
    display_options.interface_el_col = 'b';
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 1;
    anim_options.norm_val = 1e+03;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

