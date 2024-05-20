function fn_animate_result(res, steps, mod, matls, anim_options, display_options)
%FN_ANIMATE_RESULT Summary of this function goes here
%   Detailed explanation goes here

figure;

res_sum_dsps = sum(res{1,1}{1}.dsps);
anim_options.norm_val = abs(median(res_sum_dsps)); %must be positive

display_options.node_sets_to_plot(1).nd = steps{1,1}{1}.load.frc_nds;
display_options.node_sets_to_plot(2).nd = steps{1,1}{1}.mon.nds;
h_patch = fn_show_geometry(mod, matls, display_options);

fn_run_animation(h_patch, res{1,1}{1}.fld, anim_options);

end

