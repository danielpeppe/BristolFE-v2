clear;
close all;
% restoredefaultpath;
addpath("../code");

%% LOAD EXPERIMENTAL DATA

load("g4-s8_x-8000_z0_025.mat","exp_data");

%% DEFINE OPTIONS

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
%Porosity
op.porosity_plot_dists = 0;

%PARAMS
op.params = [];

%Data gen
op.data_gen = 1;
N_BATCHES = 2;
op.data_gen_batch_size = 10;
op.data_gen_load = 0;


%% RUNNING SIM

%Iterate sim for number of parameters
if op.data_gen_load
    % load('last_results.mat', 'res', 'steps', 'op_save')
    % %Plot signal
    % fn_plot_batch(op_save, res, steps, exp_data);
    % fn_plot_porosity_correlations(op_save, 'porosity')

elseif op.data_gen
    for BATCH_NUMBER = 1:N_BATCHES
        %Set default options
        [op, op_output, default_op] = fn_set_options(op);
        %Set data gen vars
        % [name, variation type, perc variation (95% of values sit val% between default op value)
        small_var = 0.01;
        med_var = 0.05;
        op.data_gen_vars = {
                    {"specimen_size", "norm", small_var}
                    {"ply0_rho_multiplier", "norm", small_var}
                    {"ply90_rho_multiplier", "norm", small_var}
                    {"ply90_D_multiplier", "norm", small_var}
                    {"ply0_D_multiplier", "norm", small_var}
                    % {"rayleigh_quality_factor", "norm", small_var} %damping changed anyway because its dependent on K and M
                    % {"interply_rho_multiplier", "norm", small_var}
                    % {"interply_D_multiplier", "norm", small_var}
                    {"water_rho_multiplier", "norm", med_var}
                    {"water_D_multiplier", "norm", med_var}
                    {"porosity", "lin", op.porosity_range}
                    {"porosity_r_sigma_tuner", "lin", [0.5 1.5]}
                    {"porosity_dist_mu_tuner", "lin", [0.5 1.5]}
                    {"porosity_dist_sigma_tuner", "lin", [0.5 1.5]}
                   };

        %Initialise variables
        n = op.data_gen_batch_size;
        res = cell(1, n);
        steps = cell(1, n);
        op_save = cell(1, n);
        %Create data folder
        dir_name = ['C:/Users/danjm/Documents/IRP/data/batch' char(string(BATCH_NUMBER))];
        if isfolder(dir_name)
            error('Directory already exists!');
        end
        mkdir(dir_name)
        
        % %Run default sim with no porosity
        % fprintf("----------------------------- No 1/%d (no porosity) ------------------------------------\n", n)
        % op.porosity = 0;
        % fn_print_default_options(op, default_op);
        % op_save{1,1} = op;
        % fprintf("---------------------------------------------------------------------------------------\n")
        % [res{1,1}, steps{1,1}] = run_sim(op, op_output);
    
        fprintf("--------------------------- RUNNNING BATCH OF SIMS ------------------------------------\n")
        for i = 1:n
            fprintf("----------------------------------- No %d/%d ------------------------------------------\n", i, n)
            op = fn_prep_data_gen(op);
            fn_print_default_options(op, default_op);
            fprintf("---------------------------------------------------------------------------------------\n")
            %Generate data
            [res{1,i}, steps{1,i}, op_save{1,i}] = run_sim(op, op_output);
            %Get attenuation
            op_save = fn_get_attenuation(op_save, res, i);
            %Save data
            dsp_sum = sum(res{1,i}{1}.dsps);
            time_data = steps{1,i}{1}.load.time;
            op_config = op_save{1,i};
            save([dir_name '/response' char(string(i)) '.mat'], 'dsp_sum', 'time_data', 'op_config')
        end
    
        %Save important plots
        fn_plot_batch(op_save, res, steps, exp_data, 1, dir_name);
        fn_plot_porosity_correlations(op_save, 'porosity', 1, dir_name)
        
        % fn_plot_porosity_correlations(op_save, 'total_n_pores', 1, dir_name)
        % fn_plot_porosity_correlations(op_save, 'porosity_r_sigma_tuner', 1, dir_name)
        % fn_plot_porosity_correlations(op_save, 'porosity_dist_mu_tuner', 1, dir_name)
        % fn_plot_porosity_correlations(op_save, 'porosity_dist_sigma_tuner', 1, dir_name)
    end
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
        % op.porosity = p;
        op.centre_freq = p;
        %%%%%%%%%%%%%%%%% Params end %%%%%%%%%%%%%%%%%
        [op, op_output] = fn_set_options(op, op_output);
        fprintf("--------------------------------------------------------------------------------------\n")
        %Get results
        [res{1,i}, steps{1,i}] = run_sim(op, op_output, {fe_options, anim_options, exp_data});
    end
    %Plot results
    fn_plot_signal(op, res, steps, exp_data)
end



