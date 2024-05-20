function [res, steps, op, mod, matls] = fn_run_sim(op, op_output, varargin)

if ~isempty(varargin)
    fe_options = varargin{1}{1};
    anim_options = varargin{1}{2};
    display_options = varargin{1}{3};
end
if op_output.animate && (fe_options.field_output_every_n_frames == inf || anim_options.repeat_n_times == 0)
    error('op_output.animate is enabled but fe_options.field_output_every_n_frames == inf, or, anim_options.repeat_n_times == 0. So animation will not be shown.')
end

%% DEFINE SPECIMEN MATERIALS

%Define rayleigh coefs
rayleigh_coefs = [0 1/(2*pi*op.centre_freq*(op.rayleigh_quality_factor * 1e4))]; %[alpha beta]

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

%Intraply boundary
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

%Define matls struct from mat
[matls, Z, R_coefs] = fn_get_matls_struct(op, mat);

%% DEFINE SHAPE OF MODEL

%Define aperture size relative to exp probe
probe_width = 0.0381;
aperture_width = double(probe_width * op.aperture_n_els/op.aperture_total_n_els);
aperture_width_8 = double(probe_width * 8/op.aperture_total_n_els);
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
el_size = fn_get_suitable_el_size(matls, op.centre_freq, op.els_per_wavelength);
%Save element size
op.el_size = el_size;
%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(model_bdry_pts, el_size);

%% PRINT IMPORTANT INFORMATION

% fprintf('el_size: %.2f um\n', el_size * 1e6)
% [max_vel, min_vel] = fn_estimate_max_min_vels(matls);
% avg_vel = (min_vel + max_vel) / 2;
% wavelengths = [min_vel max_vel avg_vel] / op.op.centre_freq * 1e6; %[um]
% fprintf('US wavelengths: min: %.2f, max: %.2f, avg: %.2f um\n', wavelengths)

%% SET MATERIALS IN MESH

%First set all elements to water
mod.el_mat_i(:) = fn_matl_i(matls,"water");
%Set specimen materials
if op.solid_specimen
    mod = fn_set_els_inside_bdry_to_mat(mod, specimen_brdy_pts, fn_matl_i(matls, op.layer1));
elseif op.composite_specimen
    if op.upper_water_present
        [mod, new_top_of_specimen] = fn_set_ply_material(mod, op, matls, specimen_brdy_pts);
    else
        %v2 does not suppot op.upper_water_present
        [mod, comp] = fn_set_ply_material_v2(mod, op, matls);
    end

    %Add porosity
    %Does not support solid specimens (because comp is needed)
    [mod, matls, op] = fn_add_porosity_v4(mod, op, matls, comp);
end

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
if strcmpi(op.src_matl,"solid_horizontal")
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
time_step = fn_get_suitable_time_step(matls, el_size, op.time_step_safety_factor);
% time_step = 4.5002e-10; %OVERRISE TIMESTEP FOR CONSISTENCY
steps{1}.load.time = 0: time_step:  op.max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, op.centre_freq, op.no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
if op.separate_receiver
    % receiver_end_pts = src_end_pts - mod.el_height*[0 1; 0 1];
    % receiver_end_pts = src_end_pts - op.specimen_size*[0 1; 0 1];
    receiver_end_pts = src_end_pts + model_width*[1 0; 1 0];
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, receiver_end_pts(1, :), receiver_end_pts(2, :), el_size / 2);
else
    steps{1}.mon.nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
end
steps{1}.mon.dfs = ones(size(steps{1}.mon.nds)) * op.src_dir;

%% DISPLAY GEOMETRY

display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
%Plot geometry
if op_output.geometry
    figure; 
    h_patch = fn_show_geometry(mod, matls, display_options);
end

%% RUN THE MODEL

if op_output.run_fea
    %Following relate to how absorbing regions are created by adding damping
    %matrix and reducing stiffness matrix to try and preserve acoustic impedance
    fe_options.damping_power_law = 3;
    fe_options.max_damping = 3.1415e+07;
    fe_options.max_stiffness_reduction = 0.01;
    %Run model
    res = fn_BristolFE_v2(mod, matls, steps, fe_options);
end

end


