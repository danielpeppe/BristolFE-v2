function op = fn_prep_data_gen(op)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(op.data_gen_vars)
    return
end

%% PREPARE VARIABLES
for i = 1:length(op.data_gen_vars)
    var_name = op.data_gen_vars{i}{1};
    var_type = op.data_gen_vars{i}{2};
    var_val = op.data_gen_vars{i}{3};
    var_default = op.(var_name);

    if strcmp(var_type, 'norm')
        mu = var_default;
        sigma = var_default*var_val;
        op.(var_name) = mu + sigma*randn;
    
    elseif strcmp(var_type, 'lin')
        up_bound = var_val(2);
        low_bound = var_val(1);
        op.(var_name) = low_bound + (up_bound - low_bound) * rand;
    end

end

