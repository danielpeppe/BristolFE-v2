function matls = fn_get_matls_struct(op, mat)
%UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
%Check composite layers are defined materials
if ~isfield(mat, op.layer1) || ~isfield(mat, op.layer2) ||...
        ~isfield(mat, op.interply_layer1) || ~isfield(mat, op.interply_layer2)
    error('Not all composite layers are defined materials')
end
%Define matls struct (just reorgonises mat into different form)
matls = struct();
matls_names = fieldnames(mat);
for i = 1:numel(matls_names)
    matls_fields = fieldnames(mat.(matls_names{i}));
    for j = 1:numel(matls_fields)
        matls(i).(matls_fields{j}) = mat.(matls_names{i}).(matls_fields{j});
    end
    matls(i).name = matls_names{i};
end
end

