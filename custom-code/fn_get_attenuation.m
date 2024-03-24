function op_save = fn_get_attenuation(op_save, res, i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Iterate over sim results
res_sum = sum(res{1,i}{1}.dsps); %tmp for readability
%Scale sim data
res_norm = res_sum/max(abs(res_sum)); %Scale exp response to match sim response
%Plot final results
n_pts_half = round(length(res_norm)/2);
op_save{1,i}.attenuation = 1/max(res_norm(n_pts_half:end));

end

