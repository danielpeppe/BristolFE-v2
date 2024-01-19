function main = fn_run_GL_sub_model(main, options)
default_options.doms_to_run = 1:numel(main.doms);
default_options.scats_to_run_in = [];
default_options.tx_to_run = 1:numel(main.steps);
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
options = fn_set_default_fields(options, default_options);
tx_to_run = options.tx_to_run;
field_output_every_n_frames = options.field_output_every_n_frames;
scats_to_run_in = options.scats_to_run_in;
doms_to_run = options.doms_to_run;


%Run main model if not already run
if ~isfield(main, 'res')
    main = fn_run_main_model(main, options);
end

%Run the scatterer models
for d = 1:numel(options.doms_to_run)
    if isempty(scats_to_run_in)
        scats_to_run = 1:numel(main.doms{doms_to_run(d)}.scats);
    else
        scats_to_run = scats_to_run_in;
    end
    for s = 1:numel(scats_to_run)
        for tx = 1:numel(tx_to_run)
            %Copy forcing from domain to scatterer model, taking account
            %the node numbers
            [   main.doms{d}.scats{s}.steps{tx}.load.frc_nds, ...
                main.doms{d}.scats{s}.steps{tx}.load.frc_dfs, ...
                main.doms{d}.scats{s}.steps{tx}.load.frcs, ...
                main.doms{d}.scats{s}.steps{tx}.load.time] = ...
            fn_GL_copy_from_main_to_local(...
                main.doms{d}.scats{s}.mod.main_nd_i, ...
                main.doms{d}.res{tx}.frc_gl_nds, ...
                main.doms{d}.res{tx}.frc_dfs, ...
                main.doms{d}.res{tx}.frcs, ...
                main.mod.time);

            %Copy monitoring from domain to scatterer model, taking account
            %the node numbers
            % main.doms{d}.scats{s}.steps{tx}.mon.nds = fn_GL_main_to_loc_nd_i(main.doms{d}.res{tx}.dsp_gl_nds, main.doms{d}.scats{s}.mod.main_nd_i);
            % main.doms{d}.scats{s}.steps{tx}.mon.dfs = main.doms{d}.res{tx}.dsp_dfs;
            % main.doms{d}.scats{s}.steps{tx}.mon.field_output_every_n_frames = options.field_output_every_n_frames;

            [   main.doms{d}.scats{s}.steps{tx}.mon.nds, ...
                main.doms{d}.scats{s}.steps{tx}.mon.dfs, ...
                main.doms{d}.scats{s}.steps{tx}.mon.field_output_every_n_frames] = ...
            fn_GL_copy_from_main_to_local(...
                main.doms{d}.scats{s}.mod.main_nd_i, ...
                main.doms{d}.res{tx}.dsp_gl_nds, ...
                main.doms{d}.scats{s}.steps{tx}.mon.dfs, ...
                options.field_output_every_n_frames);

        end
        [main.doms{d}.scats{s}.res, main.doms{d}.scats{s}.mats] = fn_BristolFE_v2(main.doms{d}.scats{s}.mod, main.matls, main.doms{d}.scats{s}.steps, options);
        %OK I think to here - correct incident field generated
        %Now need to pick up scattered field, convert to forces and
        %reciprocate!

        % for tx = 1:numel(tx_to_run)
        %     %following should be in a function (used in main model as wells)
        %     gl_gi = fn_gl_ind_for_nd_and_dof(...
        %         main.doms{d}.scats{s}.mats.gl_lookup, ...
        %         main.doms{d}.res{tx}.dsp_gl_nds, ...
        %         main.doms{d}.scats{s}.res{tx}.dsp_dfs);
        %     [main.doms{d}.scats{s}.res{tx}.frcs, force_set] = fn_convert_disps_to_forces_v2(...
        %         main.doms{d}.scats{s}.mats.K(gl_gi,gl_gi), ...
        %         main.doms{d}.scats{s}.mats.C(gl_gi,gl_gi), ...
        %         main.doms{d}.scats{s}.mats.M(gl_gi,gl_gi), ...
        %         time_step, main.doms{d}.scats{s}.res{tx}.dsps, ...
        %         main.doms{d}.res{tx}.dsp_lys, 'out');
        %     main.doms{d}.scats{s}.res{tx}.frc_nds = main.doms{d}.scats{s}.res{tx}.dsp_nds(force_set);
        %     main.doms{d}.scats{s}.res{tx}.frc_gl_nds = main.doms{d}.scats{s}.res{tx}.dsp_gl_nds(force_set);
        %     main.doms{d}.scats{s}.res{tx}.frc_dfs = main.doms{d}.scats{s}.res{tx}.dsp_dfs(force_set);        %
        % 
        % end
        % %Build global matrices for scatterer model if nesc
        % if ~isfield(main.doms{doms_to_run(d)}.scats{scats_to_run(s)}, 'mats')
        %     [main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.K, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.C, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.M, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.gl_lookup] = ...
        %     fn_build_global_matrices_v4(...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.nds, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.els, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.el_mat_i, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.el_abs_i, ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.el_typ_i, ...
        %         main.mod.matls, options);
        %     k = find(main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys2 > 0); %boundary nodes in scat model
        %     [main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_gi, ~, ~, nd_i] = ...
        %         fn_global_indices_for_all_dof_at_nodes(...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.gl_lookup, k); %assoc global DOF in scat model
        %     % main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi = main.doms{doms_to_run(d)}.mod.bdry_lys(nd_i); %%%%%
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi = ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys2(k(nd_i)); %assoc layer indices
        % end

        %Run scatterer model for each incident field - get the incident displacements at bndry nodes from main
        %**Would be more logical and tidier to put the relevant incident fields
        %into the domain models after the pristine model is run
        % bdry_gi = main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_gi;
        % tx_trans_fns = main.res.tx_rx{tx_to_run(t)}.hist(main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.main_hist_gi, :);
        % bdr_frc_gi = main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.main_hist_gi(...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi == 2 | ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi == 3);
        % time_step = main.mod.time(2) - main.mod.time(1);
        % [inc_forces, force_set] = fn_convert_disps_to_forces_v2(...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.K(bdry_gi, bdry_gi), ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.C(bdry_gi, bdry_gi), ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.M(bdry_gi, bdry_gi), ...
        %     time_step, tx_trans_fns, ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi , 'in');
        %
        % %Run the solver
        % fprintf('Running scatterer model %i in sub-doms %i for excitation %i ', s, doms_to_run(d), tx_to_run(t));
        % [scat_disp, main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.f_out, ~] = fn_explicit_dynamic_solver_v5(...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.K, main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.C, ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.M, main.mod.time, ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_gi(force_set), inc_forces, ... %forcing functions
        %     [], [], ... %disp input functions
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_gi, field_output_every_n_frames, options.use_gpu_if_available);
        %
        % %Convert outgoing displacements to forces
        % [scat_forces, force_set, ~, ~] = fn_convert_disps_to_forces_v2(...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.K(bdry_gi, bdry_gi), ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.C(bdry_gi, bdry_gi), ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mats.M(bdry_gi, bdry_gi), ...
        %     time_step, scat_disp, ...
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.bdry_lys_gi , 'out');
        % main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{t}.scat_force = scat_forces;
        %
        % %Use reciprocity to map these back to displacements at receivers
        % %NB in the final output, rx is the scattered signal + pristine result for display purposes;
        % %hist is just the scattered sgnal
        %
        % main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.hist = zeros(numel(main.res.tx_rx), numel(main.mod.time));
        % main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.rx = zeros(numel(main.res.tx_rx), numel(main.mod.time));
        %
        % %deconv_disp_input_from_main(bdr_frc_gi, :);
        % for r = 1:numel(main.res.tx_rx)
        %     rx_trans_fns = main.res.tx_rx{r}.hist(main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.mod.main_hist_gi, :);
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.hist(r,:) = ...
        %         sum(fn_convolve(...
        %         scat_forces, ...
        %         rx_trans_fns(force_set, :), 2));
        %
        %     main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.rx(r,:) = ...
        %         main.doms{doms_to_run(d)}.scats{scats_to_run(s)}.res.tx_rx{tx_to_run(t)}.hist(r,:) + ...
        %         main.res.tx_rx{t}.rx(r,:);
    end
end

end
