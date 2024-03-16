function [op, op_output] = fn_set_options(op, op_output)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Model
default_op.abs_bdry_thickness_perc = 0.2;
default_op.model_width_multiplier = 1.5;
%Resolution
default_op.els_per_wavelength = 15; %10 is default (increases are non-linear)
default_op.time_step_safety_factor = 3; %3 is default
%Signal options
default_op.centre_freq = 5e6;
default_op.no_cycles = 3;
default_op.max_time = 3.5e-6;
default_op.src_matl = 'solid'; %solid, water, solid_horizontal
default_op.src_dir = 2;
default_op.aperture_n_els = 16; %number of elements
default_op.separate_transmitter = 0; %by 1 element
default_op.separate_receiver = 0;
%Specimen
default_op.specimen_size = 4e-3; %[m]
default_op.scale_units = 1; %1000 = using mm
default_op.solid_specimen = 0;
default_op.composite_specimen = 1;
%Composite structure
default_op.ply0_orientation = 0;
default_op.ply90_orientation = 90;
default_op.n_ply_layers = 32;
default_op.n_plys_per_type = 2;
default_op.ply_symmetry = 1;
%Composite materials
default_op.layer1 = 'ply90';
default_op.layer2 = 'ply0';
%Ply options
%   density
default_op.ply0_rho_multiplier = 1;
default_op.ply90_rho_multiplier = 1;
%   stiffness
default_op.ply90_D_multiplier = 1.1;
default_op.ply0_D_multiplier = 1.1;
default_op.ply90_E_t_multiplier = 1;
default_op.ply0_E_t_multiplier = 1;
default_op.ply90_G_x_multiplier = 1;
default_op.ply0_G_x_multiplier = 1;
%   damping
default_op.rayleigh_quality_factor = 0.5; %inf disables damping
%Interply boundary options
%   v1 and v2
default_op.interply_boundary = 1;
default_op.interply_layer1 = 'resin';
default_op.interply_layer2 = 'resin';
default_op.interply_el_thickness_perc = 0.05;
%   v1
default_op.interply_midway_boundary = 1;
default_op.interply_every_layer = 1;
%   v2
default_op.interply_first_layer = 1;
default_op.interply_last_layer = 0;
default_op.interply_rho_multiplier = 1;
default_op.interply_D_multiplier = 1;

default_op.ply0b_interply_rho_multiplier = 1;
default_op.ply0b_interply_D_multiplier = 1;
default_op.ply90b_interply_rho_multiplier = 1;
default_op.ply90b_interply_D_multiplier = 1;
%   v2 - Intraply boundary options
default_op.intraply_boundary = 0;
default_op.intraply_layer1 = 'resin_intra';
default_op.intraply_layer2 = 'resin_intra';
default_op.intraply_rho_multiplier = 1;
default_op.intraply_D_multiplier = 1;
%Water boundary
default_op.water_bdry_thickness_perc = 0.25; %0.25 is default (>abs_bdry_thickness_perc)
default_op.upper_water_present = 0;
default_op.lower_water_present = 1;
default_op.water_rho_multiplier = 1.5;
default_op.water_D_multiplier = 1.5;
%Solid water options
default_op.solidwater = 0;
default_op.solidwater_rho_multiplier = 1;
default_op.solidwater_D_multiplier = 1;
%Location of transducer in water
default_op.water_interface_perc = 0; %0-1 (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
default_op.water_interface_single = 0; %0 or 1 (1 separates transducer from specimen by 1 element)
%Parameters for multiple simulations
default_op.params = [];
%Output scaling
default_op.plot_scale_dsps = 1;
default_op.plot_scale_time = 1;
%Testing horizontal speed
default_op.test_horizontal_speed = 0;
%Porosity
default_op.porosity = 0;
default_op.porosity_dist_sigma_tuner = 1.5;
default_op.porosity_dist_mu_tuner = 1; %0.5-1.5
default_op.porosity_dist_n_samples_sf = 1.5;
default_op.porosity_sigma_tuner = 1;
default_op.porosity_mu_tuner = 1;
default_op.porosity_n_pore_matls = 100;
default_op.porosity_r_min = 1e-6; %[m]
default_op.porosity_r_max = 2e-6;
default_op.porosity_r_max_multiplier = 1;
default_op.porosity_ply_matls = ["ply0","ply90"];
default_op.porosity_el_size_safety_factor = 1;
default_op.porosity_plot_dists = 0;
%Porosity Material
default_op.porosity_use_void = 0;
default_op.porosity_use_air = 0;
default_op.porosity_use_density = 0;
default_op.porosity_use_damping = 0;

%% SET DEFAULT OPTIONS

op = fn_set_default_fields(op, default_op, true); %true to print non-default options

%% VALIDATE INPUT

%turn off warning backtrace
warning('OFF', 'BACKTRACE');

%Resolution options
if op.els_per_wavelength < 15
    warning('els_per_wavelength < 15')
end
if op.time_step_safety_factor < 3
    warning('time_step_safety_factor < 3')
end
%Input options
if op.solid_specimen && op.composite_specimen
    error('Option Error: choose solid_specimen or composite_specimen')
end
if (op.lower_water_present || op.upper_water_present) && (op.abs_bdry_thickness_perc >= op.water_bdry_thickness_perc)
    warning('absorbing boundary >= water boundary')
end
%Ply options
if op.n_plys_per_type == 0
    error('Option error: n_plys_per_type must be an integer > 0')
end
if ~op.interply_boundary && op.interply_el_thickness_perc
    error('Option error: set op.interply_boundary = 0 if op.interply_el_thickness != 0')
end
%Water options
if op.water_interface_perc && op.water_interface_single
    error('Option Error: choose either op.water_interface_perc OR op.water_interface_single (not both)')
end
if (op.water_interface_perc == 1 || op.water_interface_single == 1) && op.upper_water_present == 0
    error('Option Error: op.upper_water_present = 0, set op.water_interface_perc/single = 0')
end
if (op.lower_water_present || op.upper_water_present) && ~op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc > 0')
end
if (~op.lower_water_present && ~op.upper_water_present) && op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc = 0')
end
%Solid water options
if op.solidwater && ~op.upper_water_present
    error('Option Error: set op.upper_water_present = 1 when using solidwater')
end
if (op.water_interface_perc || op.water_interface_single) && ~(strcmpi(op.src_matl,'water')) && ~op.solidwater
    error('Option Error: set op.src_matl = ''water'' when using op.water_interface_perc or op.water_interface_single')
end
%Source options
if strcmpi(op.src_matl,'water')
    op.src_dir = 4;
elseif strcmpi(op.src_matl,'solid')
    op.src_dir = 2;
elseif strcmpi(op.src_matl,'solid_horizontal')
    op.src_dir = 1;
else
    error('Option error: set op.src_matl to water, solid, or solid_horizontal')
end
%Porosity options
if  op.porosity_use_void + op.porosity_use_air + op.porosity_use_density + op.porosity_use_damping > 1
    error('Option error: chose only 1 porosity_use material')
end
if  ~op.porosity_n_pore_matls && (op.porosity_use_air || op.porosity_use_density || op.porosity_use_damping)
    error('Option error: number of pore materials is set to 0, not porosity_use materials have been set')
end
if op.porosity && ~(op.porosity_use_void || op.porosity_use_air || op.porosity_use_density || op.porosity_use_damping)
    error('Option error: porosity is enabled, but material has not been defined')
end

%% TESTING HORIZONTAL SPEED

if op.test_horizontal_speed
    op.solid_specimen = 1;
    op.src_matl = 'solid_horizontal';
    op.composite_specimen = 0;
    op.lower_water_present = 0;
    op.water_bdry_thickness_perc = 0;
    op.abs_bdry_thickness_perc = 0;
    op.model_width_multiplier = 1;
    op.rayleigh_quality_factor = inf;
    op.max_time = 2 * 3.5e-6;
end

%% VALIDATE OUTPUT

if isempty(op.params)
    op_output.plot_sim_data = 1;
    op_output.plot_exp_data = 1;
else
    op_output.justgeometry = 0;
    op_output.animate = 0;
end
if op_output.justgeometry
    fprintf("Output option enabled: JUSTGEOMETRY\n")
end
if op_output.geometry
    fprintf("Output option enabled: GEOMETRY\n")
end
if op_output.animate
    fprintf("Output option enabled: ANIMATE\n")
end

end

