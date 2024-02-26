function mod = fn_add_porosity(mod, n_pores)
%FN_GEN_VOID Summary of this function goes here
%   Detailed explanation goes here

%tmp (should be input to function)
modsize_w = max(mod.nds(:,1));
modsize_h = max(mod.nds(:,2));
modsize_avg = (modsize_w + modsize_w)/2;

%% Define porosity parameters
%Pore radius
pore_r_as_perc = 0.05;
%pore centroid locations within a region
prx_perc = 0;
pry_perc = 0.45;
prh_perc = 0.1;
prw_perc = 1;
%NOTE: can be more efficient using matrix multiplication
pore_loc_perc = rand(n_pores,2);
locx = modsize_w*(prx_perc + pore_loc_perc(:,1)*prw_perc);
locy = modsize_h*(pry_perc + pore_loc_perc(:,2)*prh_perc);


%% Select pore nodes and combine in single output
pores_out = zeros(size(mod.els,1), 1);

for i = 1:n_pores
    
    r = modsize_avg * pore_r_as_perc;
    radius = r;          % Radius of the circle
    num_points = 100;    % Number of points on the circle
    
    % Angle increments
    theta = linspace(0, 2*pi, num_points);
    
    % Cartesian coordinates (circle centered at origin)
    x = radius * cos(theta);
    y = radius * sin(theta);
    
    % Translating the circle to the new center (cx, cy)
    x_translated = x + locx(i,:);
    y_translated = y + locy(i,:);
    
    % List of points after translation
    bdry_pore = [x_translated', y_translated'];
    [in_i, out] = fn_elements_in_region(mod, bdry_pore);
    
    % combine ith pore with pore region selection
    pores_out = pores_out | in_i;
end

%label nodes with porosity as empty nodes
mod.els(pores_out, :) = [];

%Remove unused nodes
[mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);


%% plot nodes
% x = mod.nds(:,1);
% y = mod.nds(:,2);
% figure; scatter(x,y);
    
end

