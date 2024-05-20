clear;
close all;
% restoredefaultpath;
addpath("../code");

%% DEFINE OPTIONS

%Save data
op.data_gen_save_data = 0;
N_BATCHES_START = 101;
N_BATCHES_END = 101;
BATCH_SIZE = 2;

%Set data gen vars
op.porosity_range = [0 3];
%Set default options
[op, op_output, default_op] = fn_set_options(op);
% [name, variation type, perc variation (95% of values sit val% between default op value)
small_var = 0.05; med_var = 0.10; large_var = 0.20;
op.data_gen_vars = {
        {"specimen_size", "norm", small_var}
        {"ply0_rho_multiplier", "norm", med_var}
        {"ply90_rho_multiplier", "norm", med_var}
        {"ply0_D_multiplier", "norm", med_var}
        {"ply90_D_multiplier", "norm", med_var}
        {"rayleigh_quality_factor", "norm", med_var} %damping changed anyway because its dependent on K
        {"interply_rho_multiplier", "norm", med_var}
        {"interply_D_multiplier", "norm", med_var}
        {"water_rho_multiplier", "norm", large_var}
        {"water_D_multiplier", "norm", large_var}
        {"porosity", "lin", op.porosity_range}
        {"porosity_r_sigma_tuner", "lin", [0.5 1.5]}
        {"porosity_dist_mu_tuner", "lin", [0.5 1.5]}
        {"porosity_dist_sigma_tuner", "lin", [0.5 1.5]}
       };

%% RUNNING SIM

for BATCH_NUMBER = N_BATCHES_START:N_BATCHES_END
    %Initialise variables
    n = BATCH_SIZE;
    res = cell(1, n);
    steps = cell(1, n);
    op_save = cell(1, n);
    
    %Create data folder
    if op.data_gen_save_data
        batch_path = "data/batch" + BATCH_NUMBER;
        if isfolder(batch_path)
            error('Directory already exists!');
        end
        mkdir(batch_path)
    end
    
    fprintf("--------------------------- RUNNNING BATCH OF SIMS ------------------------------------\n")
    for i = 1:n
        fprintf("----------------------------------- No %d/%d ------------------------------------------\n", i, n)
        op_var = fn_prep_data_gen(op);
        fn_print_default_options(op, default_op);
        fprintf("---------------------------------------------------------------------------------------\n")

        %Generate data
        [res{1,i}, steps{1,i}, op_save{1,i}] = fn_run_sim(op_var, op_output);
        %Get attenuation
        op_save = fn_get_attenuation(op_save, res, i); %TODO: update to take op_save{1,i} as input instead of op_save

        %Save data
        if op.data_gen_save_data
            dsp_data = sum(res{1,i}{1}.dsps);
            time_data = steps{1,i}{1}.load.time;
            op_config = op_save{1,i};
            response_file_path = batch_path + "/response" + i + ".mat";
            save(response_file_path, 'dsp_data', 'time_data', 'op_config')
        end
    end

    %Save important plots
    if op.data_gen_save_data
        fn_plot_batch(op_save, res, steps, batch_path, 1);
        fn_plot_param_against_atten(op_save, 'porosity', batch_path, 1)
        % fn_plot_param_against_atten(op_save, 'total_n_pores', batch_path, 1)
        % fn_plot_param_against_atten(op_save, 'porosity_r_sigma_tuner', batch_path, 1)
        % fn_plot_param_against_atten(op_save, 'porosity_dist_mu_tuner', batch_path, 1)
        % fn_plot_param_against_atten(op_save, 'porosity_dist_sigma_tuner', batch_path, 1)
    end
end
