clear;
close all;
% restoredefaultpath;
addpath("../code");

%% LOAD EXPERIMENTAL DATA

load("g4-s8_x-8000_z0_025.mat","exp_data");

%% DEFINE OUTPUT

%Resolution options
anim_options.repeat_n_times = 10;
fe_options.field_output_every_n_frames = 50; %10 or inf is default (inf = no field output)
%Output for each sim
op_output.justgeometry = 1; %disables other outputs
op_output.geometry = 1;
op_output.run_fea = 1;
op_output.plot_sim_data = 0;
op_output.plot_exp_data = 0;
op_output.animate = 0;


%% TUNING

op.porosity = 1; %[%]
op.porosity_range = [0 5];
op.porosity_damping_tuner = 2;
%Dists
op.porosity_dist_sigma_tuner = 1;
op.porosity_dist_mu_tuner = 1; %0.5-1.5
%Matls
op.porosity_use_density = 1;
op.porosity_r_sigma_tuner = 1;
%Output
op.porosity_plot_dists = 0;

%% PARAMS

%Define input parameters
op.params = [];
op.params = [5 1 0];

%% DATA GEN

op.data_gen = 1;
op.data_gen_batch_size = 20;
op.data_gen_load = 0;
%[name, variation type, perc variation (95% of values sit val% between default op value)
small_var = 0.01;
med_var = 0.05;
large_var = 0.1;
op.data_gen_vars = {
                    % {"specimen_size", "norm", small_var}
                    % {"ply0_rho_multiplier", "norm", med_var}
                    % {"ply90_rho_multiplier", "norm", med_var}
                    % {"ply90_D_multiplier", "norm", med_var}
                    % {"ply0_D_multiplier", "norm", med_var}
                    % {"rayleigh_quality_factor", "norm", large_var}
                    % {"interply_rho_multiplier", "norm", med_var}
                    % {"interply_D_multiplier", "norm", med_var}
                    % {"water_rho_multiplier", "norm", large_var}
                    % {"water_D_multiplier", "norm", large_var}
                    {"porosity", "lin", op.porosity_range}
                    % {"porosity_r_sigma_tuner", "lin", [0.5 1.5]}
                    % {"porosity_dist_mu_tuner", "lin", [0.5 1.5]}
                    % {"porosity_dist_sigma_tuner", "lin", [0.5 1.5]}
                   };

%% RUNNING SIM

%Iterate sim for number of parameters
if op.data_gen_load
    load('last_results.mat', 'res', 'steps', 'op_save')
    %Plot signal
    op_save = fn_check_porosity_data(op_save, res, steps, exp_data);
elseif op.data_gen
    %Set default options
    [op, op_output, default_op] = fn_set_options(op);

    %Initialise variables
    n = op.data_gen_batch_size + 1;
    res = cell(1, n);
    steps = cell(1, n);
    op_save = cell(1, n);
    
    %Run default sim with no porosity
    fprintf("----------------------------- No 1/%d (no porosity) ------------------------------------\n", n)
    op.porosity = 0;
    fn_print_default_options(op, default_op);
    op_save{1,1} = op;
    fprintf("---------------------------------------------------------------------------------------\n")
    [res{1,1}, steps{1,1}] = run_sim(op, op_output);

    %Run Batch
    fprintf("--------------------------- RUNNNING BATCH OF SIMS ------------------------------------\n")
    for i = 2:n
        fprintf("----------------------------------- No %d/%d ------------------------------------------\n", i, n)
        op = fn_prep_data_gen(op);
        fn_print_default_options(op, default_op);
        op_save{1,i} = op;
        fprintf("---------------------------------------------------------------------------------------\n")
        %Generate data
        [res{1,i}, steps{1,i}] = run_sim(op, op_output);
        
        %Save data
        %if op.data_check, then save geometries
        % figure;
        % plot(steps{1,i}{1}.load.time, sum(res{1,i}{1}.dsps));
        % xlabel("Time (s)")
        % ylabel("Magnitude (-)")
    
        %TODO: Check Batch
    end

    %Plot signal
    op_save = fn_check_porosity_data(op_save, res, steps, exp_data);

    %Save data
    save('last_results.mat', 'res', 'steps', 'op_save')

elseif isempty(op.params)
    fprintf("--------------------------- RUNNNING ONE SIM -----------------------------------\n")
    [op, op_output] = fn_set_options(op, op_output);
    fprintf("--------------------------------------------------------------------------------\n")
    %Get results
    [res{1}, steps{1}] = run_sim(op, op_output, {fe_options, anim_options, exp_data});

else
    fprintf("--------------------------- RUNNNING MULTIPLE SIMS -----------------------------------\n")
    n = length(op.params);
    res = cell(1, n);
    steps = cell(1, n);
    for i = 1:n
        fprintf("----------------------------------- No %d/%d ------------------------------------------\n", i, n)
        p = op.params(i);
        %%%%%%%%%%%%%%%% Params start %%%%%%%%%%%%%%%%
        op.porosity = p;
        %%%%%%%%%%%%%%%%% Params end %%%%%%%%%%%%%%%%%
        [op, op_output] = fn_set_options(op, op_output);
        fprintf("--------------------------------------------------------------------------------------\n")
        %Get results
        [res{1,i}, steps{1,i}] = run_sim(op, op_output, {fe_options, anim_options, exp_data});
    end
    %Plot results
    fn_plot_signal(op, res, steps, exp_data)
end

p = zeros(size(op_save));
a = zeros(size(op_save));
for i = 1:numel(op_save)
    p(i) = op_save{i}.porosity;
    a(i) = op_save{i}.attenuation;
end
[p_sort, p_i] = sort(p);
a_sort = a(p_i);

figure; plot(p_sort, a_sort)
xlabel('Porosity [%]')
ylabel('Attenuation [% back-wall echo]')

