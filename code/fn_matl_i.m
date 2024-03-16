function matl_i = fn_matl_i(matls, matl_string)
%FN_MATL_I Summary of this function goes here
%   Detailed explanation goes here
matl_i = [];
if ~(isstring(matl_string) || ischar(matl_string))
    error('Input must be name of material as a string')
end
for i = 1:numel(matls)
    if strcmpi(matls(i).name, matl_string)
        matl_i = i;
    end
end
if isempty(matl_i)
    error('Material index could not be found')
end
end

