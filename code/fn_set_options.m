function op = fn_set_options(op)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Model
default_op.abs_bdry_thickness_perc = 0.2; %0.2 is default (relative to specimen_size)
default_op.model_size_w_multiplier = 2.52; %relative to aperture size
%Resolution
default_op.els_per_wavelength = 30; %10 is default (increases are non-linear)
default_op.time_step_safety_factor = 3; %3 is default
%Signal options
default_op.centre_freq = 5e6;
default_op.horizontal_src = 0;
default_op.aperture_n_els = 8; %number of elements
default_op.separate_transmitter = 0; %by 1 element
default_op.separate_receiver = 0;
%Specimen
default_op.specimen_size = 4e-3; %[m]
default_op.scale_model = 1;
default_op.solid_specimen = 0;
default_op.composite_specimen = 1;
%Composite structure
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
default_op.ply90_D_multiplier = 1;
default_op.ply0_D_multiplier = 1;
default_op.ply90_shear_wave_velocity_by_E_t = 1; %1 is default
default_op.ply0_shear_wave_velocity_by_E_t = 1; %1 is default
%   damping
default_op.rayleigh_quality_factor = inf; %inf disables damping
%Interply boundary options
%   v1 and v2
default_op.interply_layer1 = 'resin';
default_op.interply_layer2 = 'resin';
default_op.interply_boundary = 1; %1 is default
%   v1
default_op.interply_midway_boundary = 1;
default_op.interply_every_layer = 1;
%   v2
default_op.interply_first_layer = 1;
default_op.interply_last_layer = 1;
default_op.interply_rho_multiplier = 1;
default_op.interply_D_multiplier = 1;
%   v2 - Intraply boundary options
default_op.intraply_layer1 = 'resinb';
default_op.intraply_layer2 = 'resinb';
default_op.intraply_rho_multiplier = 1;
default_op.intraply_D_multiplier = 1;
%Water boundary
default_op.water_bdry_thickness_perc = 0.25; %0.25 is default (>abs_bdry_thickness_perc)
default_op.upper_water_present = 0;
default_op.lower_water_present = 1;
default_op.water_rho_multiplier = 1;
default_op.water_D_multiplier = 1;
%Solid water options
default_op.solidwater = 0;
default_op.solidwater_rho_multiplier = 1;
default_op.solidwater_D_multiplier = 1;
%Location of transducer in water
default_op.water_interface_perc = 0; %0-1 (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
default_op.water_interface_single = 0; %0 or 1 (1 separates transducer from specimen by 1 element)

%% Set default options
op = fn_set_default_fields(op, default_op, true); %true to print non-default options

%% Validate input

%Resolution options
if op.els_per_wavelength < 30
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('CAUTION: els_per_wavelength < 30')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
elseif op.time_step_safety_factor < 3
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('CAUTION: time_step_safety_factor < 3')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%Input options
if op.solid_specimen && op.composite_specimen
    error('Option Error: choose solid_specimen or composite_specimen')
elseif (op.lower_water_present || op.upper_water_present) && (op.abs_bdry_thickness_perc >= op.water_bdry_thickness_perc)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('CAUTION: absorbing boundary >= water boundary')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
%Ply options
if op.n_plys_per_type == 0
    error('Option error: n_plys_per_type must be an integer > 0')
end
%Water options
if op.water_interface_perc && op.water_interface_single
    error('Option Error: choose either op.water_interface_perc OR op.water_interface_single (not both)')
elseif (op.water_interface_perc == 1 || op.water_interface_single == 1) && op.upper_water_present == 0
    error('Option Error: op.upper_water_present == 0, set op.water_interface_perc/single == 0')
elseif (op.lower_water_present || op.upper_water_present) && ~op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc > 0')
elseif (~op.lower_water_present && ~op.upper_water_present) && op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc = 0')
end
%Solid water options
if op.solidwater && ~op.upper_water_present
    error('Option Error: set op.upper_water_present == 1 when using solidwater')
end
end

