function [matls, Z, R_coefs] = fn_get_matls_struct(op, mat, varargin)
%UNTITLED Summary of this function goes here
% %   Detailed explanation goes here

%Check if first time running
if isempty(varargin)
    matls = struct();
    n_matls_old = 0;
    %Check composite layers are defined materials
    if ~isfield(mat, op.layer1) || ~isfield(mat, op.layer2) ||...
            ~isfield(mat, op.interply_layer1) || ~isfield(mat, op.interply_layer2)
        error('Not all composite layers are defined materials')
    end
else
    matls = varargin{1};
    matls_names_old = {matls(:).name}';
    n_matls_old = size(matls_names_old, 1);
end

matls_names = fieldnames(mat);
n_matls = numel(matls_names);
for i = 1:n_matls
    %Reorganise into matls struct
    matls_fields = fieldnames(mat.(matls_names{i}));
    for j = 1:numel(matls_fields)
        matls(i + n_matls_old).(matls_fields{j}) = mat.(matls_names{i}).(matls_fields{j});
    end
    matls(i + n_matls_old).name = matls_names{i};
    %Check fields are named correctly
    if ~isfield(matls,'rho') || ~isfield(matls,'D') || ~isfield(matls,'name') ||...
        ~isfield(matls,'rayleigh_coefs') || ~isfield(matls,'col') || ~isfield(matls,'el_typ')
        warning('OFF', 'BACKTRACE');
        warning('Unknown material field (not rho, D, rayleigh_coefs, col, or el_typ)')
    end
        
    % %Scale values according to scale of model
    % if isfield(matls(i), 'rho')
    %     matls(i).rho = matls(i).rho / op.scale_units^3;
    % end
    % if isfield(matls(i), 'D')
    %     matls(i).D = matls(i).D / op.scale_units^2;
    % end
    % if isfield(matls(i), 'rayleigh_coefs')
    %     matls(i).rayleigh_coefs = matls(i).rayleigh_coefs / op.scale_units;
    % end
end


%Return accoustic impedances
Z = zeros(n_matls, 1);
for i = 1:n_matls
    [max_vel_i, min_vel_i] = fn_estimate_max_min_vels(matls(i), op.scale_units);
    avg_vel = (max_vel_i + min_vel_i)/2;
    Z(i) = matls(i).rho * avg_vel;
end
%Calculate reflectivity coefficients (from i into j)
R_coefs = zeros(n_matls, n_matls);
for i = 1:n_matls
    for j = 1:n_matls
        R_coefs(i, j) = (Z(j) - Z(i))/(Z(j) + Z(i));
    end
end

end
