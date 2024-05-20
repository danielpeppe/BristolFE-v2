function op_save = fn_get_attenuation(op_save, res, i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Iterate over sim results
res_sum = sum(res{1,i}{1}.dsps); %tmp for readability
%Plot final results
n_pts_half = round(length(res_sum)/2);
amplitude_drop = max(abs(res_sum))/max(res_sum(n_pts_half:end));
op_save{1,i}.attenuation = 20*log10(amplitude_drop) / (2*op_save{i}.specimen_size*1000);

end

