function matl_i_arr = fn_matl_i(matls, matls_name_arr)
%FN_MATL_I Summary of this function goes here
%   Detailed explanation goes here
if ischar(matls_name_arr)
    error('Name of materials has to be a string')
end

[~, i] = ismember({matls(:).name}, matls_name_arr);
matl_i_arr = find(i > 0);
if isempty(matl_i_arr)
    error('Material index could not be found')
end
end

