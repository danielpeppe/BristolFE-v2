function op = fn_set_options(op)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Model
default_op.abs_bdry_thickness_perc = 0.2; %0.2 is default (relative to specimen_size)
default_op.model_size_w_multiplier = 1.5;
%Location of transducer in water
default_op.upper_water_present = 0;
default_op.lower_water_present = 1;
%Water boundary
default_op.water_bdry_thickness_perc = 0.2; %0.2 is default (>abs_bdry_thickness_perc)
%Water interface
default_op.water_interface_perc = 0; %0-1 (1 separates transducer from specimen by water_boundary_thickness) (if you want src in material, set to 0 and manually edit src_offset)
default_op.water_interface_single = 0; %0 or 1 (1 separates transducer from specimen by 1 element)
%Specimen
default_op.specimen_size = 4e-3; %[mm]
default_op.solid_specimen = 0;
default_op.composite_specimen = 1;
%Composite structure
default_op.n_ply_layers = 32;
default_op.n_plys_per_type = 2;
default_op.ply_symmetry = 1;
%Ply options
default_op.wave_velocity_by_E_t = 1; %1 is default (adjusts E_t stiffness)
default_op.back_wall_reflection_by_water_density = 1; %1 is default
default_op.boundary_density_multiplier = 1;
default_op.boundary_stiffness_multiplier = 1;
%Signal options
default_op.horizontal_src = 0;
%Output
default_op.geometry = 0;
default_op.run_fea = 1;
default_op.plot_sim_data = 1;
default_op.plot_exp_data = 1;
default_op.animate = 1;

op = fn_set_default_fields(op, default_op);

%% Validate input
%Input options
if op.solid_specimen && op.composite_specimen
    error('Option Error: choose solid_specimen or composite_specimen')
end
%Water options
if op.water_interface_perc && op.water_interface_single
    error('Option Error: Choose either op.water_interface_perc OR op.water_interface_single (not both)')
elseif xor(op.water_interface_perc || op.water_interface_single, op.upper_water_present)
    error('Option Error: set op.water_interface_perc/single correctly')
elseif (op.lower_water_present || op.upper_water_present) && ~op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc > 0')
elseif (~op.lower_water_present && ~op.upper_water_present) && op.water_bdry_thickness_perc
    error('Option Error: set op.water_bdry_thickness_perc = 0')
end
%Output options
if (op.animate || op.plot_sim_data) && ~op.run_fea
    error('Option Error: set run_fea = 1 for results')
end

%% If Geometry = 1 then just plot geometry
if op.geometry
    op.run_fea = 0;
    op.plot_sim_data = 0;
    op.plot_exp_data = 0;
    op.animate = 0;
end

end

