function [op, op_output, default_op] = fn_set_options(op, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    op_output = varargin{1};
else
    op_output.geometry = 0;
    op_output.run_fea = 1;
    op_output.plot_sim_data = 0;
    op_output.plot_exp_data = 0;
    op_output.animate = 0;
end

%Experimental data
default_op.load_exp = 0;
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
default_op.src_matl = "solid"; %solid, water, solid_horizontal
default_op.src_dir = 2;
default_op.aperture_total_n_els = 128;
default_op.aperture_n_els = 16; %number of elements
default_op.separate_transmitter = 0; %by 1 element
default_op.separate_receiver = 0;
%Specimen
default_op.specimen_size = 4e-3; %[m]
default_op.solid_specimen = 0;
default_op.composite_specimen = 1;
%Composite structure
default_op.ply0_orientation = 0;
default_op.ply90_orientation = 90;
default_op.n_ply_layers = 32;
default_op.n_plys_per_type = 2;
default_op.ply_symmetry = 1;
%Composite materials
default_op.layer1 = "ply90";
default_op.layer2 = "ply0";
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
default_op.interply_layer1 = "resin";
default_op.interply_layer2 = "resin";
default_op.interply_el_thickness_perc = 0.05;
%   v1
default_op.interply_midway_boundary = 1;
default_op.interply_every_layer = 1;
%   v2
default_op.interply_first_layer = 1;
default_op.interply_last_layer = 0;
default_op.interply_rho_multiplier = 1;
default_op.interply_D_multiplier = 1;
%   v2 - Intraply boundary options
default_op.intraply_boundary = 0;
default_op.intraply_layer1 = "resin_intra";
default_op.intraply_layer2 = "resin_intra";
default_op.intraply_rho_multiplier = 1;
default_op.intraply_D_multiplier = 1;
%Water boundary
default_op.water_bdry_thickness_perc = 0.25; %0.25 is default (>abs_bdry_thickness_perc)
default_op.upper_water_present = 0;
default_op.lower_water_present = 1;
default_op.water_rho_multiplier = 1.5;
default_op.water_D_multiplier = 1.5;
%Location of transducer in water
default_op.water_interface_perc = 0; %0-1 (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
default_op.water_interface_single = 0; %0 or 1 (1 separates transducer from specimen by 1 element)
%Parameters for multiple simulations
default_op.params = [];
%Testing horizontal speed
default_op.test_horizontal_speed = 0;
%Porosity
default_op.porosity = 0;
default_op.porosity_matls_for_pores = ["ply0" "ply90"]; %only ply materials supported (due to n_pores calculation)
%   Spatial
default_op.porosity_dist_sigma_tuner = 1;
default_op.porosity_dist_mu_tuner = 1; %0.5-1.5
default_op.porosity_dist_n_samples_sf = 3;
%   Radii
default_op.porosity_r_sigma_tuner = 1;
default_op.porosity_n_pore_matls = 100;
default_op.porosity_r_min = 1e-6; %[m]
default_op.porosity_r_max = 5e-6;
%   Voids
default_op.porosity_use_void = 0;
%   Density
default_op.porosity_use_density = 1;
%   Damping
default_op.porosity_damping_tuner = 3.2;
%   Output
default_op.porosity_plot_dists = 0;
%Data generation options
default_op.data_gen_vars = {};
default_op.data_gen_save_data = 0;
default_op.porosity_range = [0 3];

%% SET DEFAULT OPTIONS

op = fn_set_default_fields(op, default_op);
%Print which options are changed from default
fn_print_default_options(op, default_op)

%% VALIDATE INPUT

%turn off warning backtrace
warning("OFF", "BACKTRACE");

%Experimental data options
if ~op.load_exp && op_output.plot_exp_data
    error("op.load_exp is not enabled so op_output.plot_exp_data is not possible")
end
%Resolution options
if op.els_per_wavelength ~= 15
    warning("els_per_wavelength may be overrided in run_sim, please check!")
elseif op.els_per_wavelength < 15
    warning("els_per_wavelength < 15")
end
if op.time_step_safety_factor < 3
    warning("time_step_safety_factor < 3")
end
%Input options
if op.solid_specimen && op.composite_specimen
    error("Option Error: choose solid_specimen or composite_specimen")
end
if (op.lower_water_present || op.upper_water_present) && (op.abs_bdry_thickness_perc >= op.water_bdry_thickness_perc)
    warning("absorbing boundary >= water boundary")
end
%Ply options
if op.n_plys_per_type == 0
    error("Option error: n_plys_per_type must be an integer > 0")
end
if ~op.interply_boundary && op.interply_el_thickness_perc
    error("Option error: set op.interply_boundary = 0 if op.interply_el_thickness_perc != 0")
end
%Water options
if op.water_interface_perc && op.water_interface_single
    error("Option Error: choose either op.water_interface_perc OR op.water_interface_single (not both)")
end
if (op.water_interface_perc == 1 || op.water_interface_single == 1) && op.upper_water_present == 0
    error("Option Error: op.upper_water_present = 0, set op.water_interface_perc/single = 0")
end
if (op.lower_water_present || op.upper_water_present) && ~op.water_bdry_thickness_perc
    error("Option Error: set op.water_bdry_thickness_perc > 0")
end
if (~op.lower_water_present && ~op.upper_water_present) && op.water_bdry_thickness_perc
    error("Option Error: set op.water_bdry_thickness_perc = 0")
end
%Source options
if strcmpi(op.src_matl, "water")
    op.src_dir = 4;
elseif strcmpi(op.src_matl, "solid")
    op.src_dir = 2;
elseif strcmpi(op.src_matl, "solid_horizontal")
    op.src_dir = 1;
else
    error("Option error: set op.src_matl to water, solid, or solid_horizontal")
end
%Porosity options
if ischar(op.porosity_matls_for_pores)
    error("Option error: op.porosity_ply_matls has to be a list of strings (not chars)")
end
if  op.porosity_use_void + op.porosity_use_density > 1
    error("Option error: chose only 1 porosity_use material")
end
if  ~(op.porosity_n_pore_matls || op.porosity) && (op.porosity_use_void || op.porosity_use_density)
    error("Option error: number of pore materials is set to 0, not porosity_use materials have been set")
end
if op.porosity == 1 && ~(op.porosity_use_void || op.porosity_use_density)
    error("Option error: porosity material needs to be defined or porosity needs to be turned off")
end

%% TESTING HORIZONTAL SPEED

if op.test_horizontal_speed
    op.solid_specimen = 1;
    op.src_matl = "solid_horizontal";
    op.src_dir = 1;
    op.composite_specimen = 0;
    op.lower_water_present = 0;
    op.water_bdry_thickness_perc = 0;
    op.abs_bdry_thickness_perc = 0;
    op.model_width_multiplier = 1;
    op.rayleigh_quality_factor = inf;
    op.max_time = 2 * 3.5e-6;
end

%% OUTPUT OPTIONS

if op_output.geometry
    fprintf("Output option enabled: GEOMETRY\n")
end
if op_output.animate
    fprintf("Output option enabled: ANIMATE\n")
end

end

