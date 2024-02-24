clear;
close all;
restoredefaultpath;
addpath('../code');

%% KEY MODELLING PARAMETERS

%%%%% Current Working Model
%   1. upper_water_present = 0 - this ensures the instability at water
%   boundary does not occure at transducer. Therefore, instablity
%   'reaches' src later on in simulation

%Resolution
els_per_wavelength = 2; %10 is default (increases are non-linear)
time_step_safety_factor = 3; %3 is default

%Specimen size and absorbing boundary
specimen_size = 4e-3; %[mm]
abs_bdry_thickness_perc = 0.1; %0.1 is default (relative to specimen_size)
%Material indices (Ply Material layers = 1 and 2) (indices cannot be skipped)
ply90_matl_i = 1;
ply0_matl_i = 2;
water_matl_i = 3;
steel_matl_i = 4; %DEBUGGING

%Location of transducer in water
op.upper_water_present = 0;
op.lower_water_present = 1;
%Water boundary
water_bdry_thickness_perc = 0.2; %0.2 is default (>abs_bdry_thickness_perc)
%Water interface
op.water_interface_perc = 0; %0-1 (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
op.water_interface_single = 0; %0 or 1 (1 separates transducer from specimen by 1 element)

%Setup
op.steel_and_water = 0;
op.plys_and_water = 1;
%Plys
op.n_ply_layers = 2;
op.n_plys_per_type = 1;
op.ply_symmetry = 0;
%Signal
max_time = 5e-6;
%FEA options
fe_options.field_output_every_n_frames = 100; %10 or inf is default (inf = no field output)
fe_options.use_gpu_if_present = 1;
%Output
op.geometry = 1;
if ~op.geometry %tmp
    op.run_fea = 1;
    op.plot_sim_data = 1;
    op.plot_exp_data = 1;
    op.animate = 1;
else
    op.run_fea = 0;
    op.plot_sim_data = 0;
    op.plot_exp_data = 0;
    op.animate = 0;
end
%% CALCULATE MORE PARAMETERS (CHECK OPTIONS)

%Water options
if op.water_interface_perc && op.water_interface_single
    error('Option Error: Choose either op.water_interface_perc OR op.water_interface_single (not both)')
elseif xor(op.water_interface_perc || op.water_interface_single, op.upper_water_present)
    error('Option Error: set op.water_interface_perc/single correctly')
elseif (op.lower_water_present || op.upper_water_present) && ~water_bdry_thickness_perc
    error('Option Error: set water_bdry_thickness_perc > 0')
elseif (~op.lower_water_present && ~op.upper_water_present) && water_bdry_thickness_perc
    error('Option Error: set water_bdry_thickness_perc = 0')
end
%Output options
if (op.animate || op.plot_sim_data || op.plot_exp_data) && ~op.run_fea
    error('Option Error: set run_fea = 1 for results')
end

%% DEFINE MATERIAL

%Steel (FOR DEBUGGING)
matls(steel_matl_i).rho = 8900;
matls(steel_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
matls(steel_matl_i).col = hsv2rgb([3/4,0.5,0.80]);
matls(steel_matl_i).name = 'Steel';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Ply Material properties
E_fib = 163e9; %[GPa]
G_fib = 5e9;
v_fib = 0.32;
E_t = 12e9;
G_t = 3.98e9;
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
model_size_h = specimen_size + water_brdy_thickness*(op.upper_water_present + op.lower_water_present);
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
    0,            wbt*op.lower_water_present
    model_size_w, wbt*op.lower_water_present
    model_size_w, wbt*op.lower_water_present + specimen_size
    0,            wbt*op.lower_water_present + specimen_size];
%Define top of specimen for later use
top_of_specimen = specimen_brdy_pts(3,2);

%Define start of absorbing boundary region and its thickness
abt = abs_bdry_thickness; %tmp for readability
abs_bdry_pts = [
    abt,                abt*op.lower_water_present
    model_size_w - abt, abt*op.lower_water_present
    model_size_w - abt, model_size_h - op.upper_water_present*abt
    abt,                model_size_h - op.upper_water_present*abt];

%% DEFINE SIGNAL

centre_freq = 5e6;
no_cycles = 4;
if op.upper_water_present
    src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)
else
    src_dir = 2;
end

%% DEFINE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%% DEFINE MATERIALS

%First set all elements to water
mod.el_mat_i(:) = water_matl_i;
%Set specimen materials
if op.steel_and_water
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, steel_matl_i);
end
if op.plys_and_water
    %Set ply materials (input 1 = layer 1, input 2 = layer 2)
    [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, 1, 2, specimen_brdy_pts, model_size_h);
end
%Add interface elements
mod = fn_add_fluid_solid_interface_els(mod, matls);
%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%% ADD POROSITY (WIP)
% n_pores = 10;
% mod = fn_add_porosity(mod, n_pores);

%% DEFINE TRANSDUCER

%Define a line along which sources will be placed to excite waves
%Adjust src offset if upper water is present
if op.upper_water_present
    %Redefine top of specimen
    top_of_specimen = new_top_of_specimen;
    if op.water_interface_single
        src_offset = mod.el_height;
    elseif op.water_interface_perc
        src_offset = op.water_interface_perc*wbt;
    else
        error('Critical Option Error: problem if water options')
    end
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

if op.geometry
    figure; 
    display_options.interface_el_col = 'b';
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
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
    fprintf('max dsp: %e, average dsp: %e median dsp: %e\n',max(res_sum_dsps), mean(res_sum_dsps), median(res_sum_dsps))
end

%% SHOW THE RESULTS

if op.plot_sim_data
    %Show the history output as a function of time - here we just sum over all 
    %the nodes where displacments were recorded
    figure;
    plot(steps{1}.load.time, res_sum_dsps);
    xlabel('Time (s)')
    hold on
    %Plot experimental data on top
    if op.plot_exp_data
        load('g4-s8_x-8000_z0_025.mat');
        aperture = 1:8;
        i = ismember(exp_data.tx, aperture) & ismember(exp_data.rx,aperture);
        time_data_i = sum(exp_data.time_data(:,i),2);
        scale_dsp = max(res_sum_dsps)/max(time_data_i); %Scale exp response to match sim response
        translate_time = 1.194e-05; %Start of exp response
        plot(exp_data.time - translate_time, scale_dsp*time_data_i);
    end
    hold off
end
%op.animate result
if op.animate
    figure;
    display_options.interface_el_col = 'b';
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 10;
    anim_options.norm_val = median(res_sum_dsps);
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

