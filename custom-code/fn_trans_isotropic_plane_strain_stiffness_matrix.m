function stiffness_matrix = fn_trans_isotropic_plane_strain_stiffness_matrix(ply_orientation,E_fib,G_fib,v_fib,E_t,G_t,v_t)
%SUMMARY
%   Returns 3 x 3 plane strain stiffness matrix for transversely isotropic material
%   relating stress [sigma_xx, sigma_yy, sigma_xy] to strain [epsilon_xx,
%   epsilon_yy, epsilon_xy].
%INPUTS
%   orientation - orientation of composite ply layer (0 or 90)
%   E_fib - Young's modulus in fibre direction5e+9
%   G_fib - Shear modulus in 'fibre' direction
%   v_fib - Poisson's ratio perpendicular to fibre plane
%   E_t - Young's modulus in transverse direction
%   G_t - Shear modulus in 'transverse' direction
%   v_t - Poisson's ratio perpendiular to transverse plane
%OUTPUTS
%   stiffness_matrix - 3 x 3 stiffness matrix
%--------------------------------------------------------------------------

%input error checks
if ~isscalar(ply_orientation) || ~isscalar(E_fib) || ~isscalar(G_fib) ||...
        ~isscalar(v_fib) || ~isscalar(E_t) || ~isscalar(G_t) || ~isscalar(v_t)
    error('input to fn_trans_isotropic_plane_strain_stiffness_matrix has to be scalar');
end
if ~(ply_orientation == 0 || ply_orientation == 90)
    error("orientations not equal to 0 or 90 are not supported")
end

if ply_orientation == 0
    % start with 6x6 compliance matrix
    compliance_matrix = [1/E_t,         -v_t/E_t,       -v_fib/E_fib, 0, 0, 0;...
                         -v_t/E_t,      1/E_t,          -v_fib/E_fib, 0, 0, 0;...
                         -v_fib/E_fib,  -v_fib/E_fib,   1/E_fib,      0, 0, 0;...
                         0,             0,              0,            1/(2*G_fib), 0, 0;...
                         0,             0,              0,            0, 1/(2*G_fib), 0;...
                         0,             0,              0,            0, 0, 1/(2*G_t)];

elseif ply_orientation == 90
    compliance_matrix = [1/E_fib,       -v_fib/E_fib,   -v_fib/E_fib, 0, 0, 0;...
                         -v_fib/E_fib,  1/E_t,          -v_t/E_t,     0, 0, 0;...
                         -v_fib/E_fib,  -v_t/E_t,       1/E_t,        0, 0, 0;...
                         0,             0,              0,            1/(2*G_t), 0, 0;...
                         0,             0,              0,            0, 1/(2*G_fib), 0;...
                         0,             0,              0,            0, 0, 1/(2*G_fib)];
else
    error("ply orientations not equal to 0 or 90 are not supported")
end

% invert compliance matrix to get stiffness matrix
stiffness_matrix = inv(compliance_matrix);
% remove z-axis component for plane strain
stiffness_matrix(3:5,:) = [];
stiffness_matrix(:,3:5) = [];

end

