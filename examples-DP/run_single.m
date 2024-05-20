clear;
close all;
% restoredefaultpath;
addpath("../code");

%% DEFINE OPTIONS

%Resolution options
anim_options.repeat_n_times = 10;
fe_options.field_output_every_n_frames = 50; %10 or inf is default (inf = no field output)
%Output for each sim
op_output.geometry = 0;
op_output.run_fea = 1;
op_output.plot_sim_data = 0;
op_output.plot_exp_data = 0;
op_output.animate = 0;

%Geometry display options
display_options.interface_el_col = hsv2rgb([0.50,.75,.60]);
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).col = 'b.';

%% LOAD EXPERIMENTAL DATA (OPTIONAL)

op.load_exp = 0;
if op.load_exp
    %Load exp data
    load("g4-s8_x-8000_z0_025.mat","exp_data");
else
    exp_data = [];
end

%% RUNNING SIM

fprintf("--------------------------- RUNNNING ONE SIM -----------------------------------\n")
[op, op_output] = fn_set_options(op, op_output);
fprintf("--------------------------------------------------------------------------------\n")
op_save = {}; %Initialise op cell to store final options configured for sim after run
%Get results
[res{1,1}, steps{1,1}, op_save{1,1}, mod, matls] = fn_run_sim(op, op_output, {fe_options, anim_options, display_options});
%Plot signal
if op_output.plot_sim_data || op_output.plot_exp_data
    fn_plot_results(op, op_output, res, steps, exp_data)
end
%Plot result (only works for single run)
if op_output.animate
    fn_animate_result(res, steps, mod, matls, anim_options, display_options)
end



