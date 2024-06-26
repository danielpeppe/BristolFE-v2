function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE3_v2(nds, els, D, rho, rayleigh_coefs, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_CPE3_v2.m.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[tdofs_to_use = [] - optional string listing the DoFs to use, e.g. '12'. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (11-Jan-2024 16:00:04)

%Deal with optional argument about which DOFs to use
if isempty(varargin)
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end

%Record the local node numbers of the element stiffness matrices
loc_nd = [1  1  1  2  2  2  3  3  3];

%Record the local DOFs of the element stiffness matrices
loc_df = [1  2  3  1  2  3  1  2  3];

%Get the DOFs if not specified
if isempty(dofs_to_use)
	dofs_to_use = unique(loc_df);
end

%If any inputs blank, return at this point with just the loc_nd and loc_df
if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)
	el_K = [];
	el_M = [];
    el_C = [];
	el_Q = [];
	[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);
	return
end

%Temporary matrices of nodal coordinates to save time
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);
nds_3_1 = nds(els(:, 3), 1);
nds_3_2 = nds(els(:, 3), 2);

%Jacobian
J = zeros(size(els, 1), 1, 1);
J(:, 1, 1) = nds_1_1 .* nds_2_2 - nds_1_2 .* nds_2_1 - nds_1_1 .* nds_3_2 + nds_1_2 .* nds_3_1 + nds_2_1 .* nds_3_2 - nds_2_2 .* nds_3_1;

%Stiffness matrix
el_K = zeros(size(els, 1), 9, 9);
el_K(:, 1, 1) = (((D(1, 1) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 1) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(1, 6) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 1, 2) = (((D(1, 6) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(1, 2) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 2) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 1, 3) = (((D(1, 5) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 5) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(1, 4) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 4) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 1, 4) = (((D(1, 6) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(1, 1) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 1) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 1, 5) = (((D(1, 2) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 2) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(1, 6) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 1, 6) = (((D(1, 4) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 4) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(1, 5) .* (nds_2_2 - nds_3_2)) ./ J - (D(6, 5) .* (nds_2_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 1, 7) = ((nds_1_2 - nds_2_2) .* (D(1, 1) .* nds_2_2 - D(1, 1) .* nds_3_2 - D(6, 1) .* nds_2_1 + D(6, 1) .* nds_3_1)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(1, 6) .* nds_2_2 - D(1, 6) .* nds_3_2 - D(6, 6) .* nds_2_1 + D(6, 6) .* nds_3_1)) ./ (2 .* J);
el_K(:, 1, 8) = ((nds_1_2 - nds_2_2) .* (D(1, 6) .* nds_2_2 - D(1, 6) .* nds_3_2 - D(6, 6) .* nds_2_1 + D(6, 6) .* nds_3_1)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(1, 2) .* nds_2_2 - D(1, 2) .* nds_3_2 - D(6, 2) .* nds_2_1 + D(6, 2) .* nds_3_1)) ./ (2 .* J);
el_K(:, 1, 9) = ((nds_1_2 - nds_2_2) .* (D(1, 5) .* nds_2_2 - D(1, 5) .* nds_3_2 - D(6, 5) .* nds_2_1 + D(6, 5) .* nds_3_1)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(1, 4) .* nds_2_2 - D(1, 4) .* nds_3_2 - D(6, 4) .* nds_2_1 + D(6, 4) .* nds_3_1)) ./ (2 .* J);
el_K(:, 2, 1) = (((D(2, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(2, 1) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 1) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 2, 2) = (((D(2, 2) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 2) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(2, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 2, 3) = (((D(2, 4) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 4) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(2, 5) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 5) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 2, 4) = (((D(2, 1) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 1) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(2, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 2, 5) = (((D(2, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(2, 2) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 2) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 2, 6) = (((D(2, 5) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 5) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(2, 4) .* (nds_2_1 - nds_3_1)) ./ J - (D(6, 4) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 2, 7) = ((nds_1_1 - nds_2_1) .* (D(2, 6) .* nds_2_1 - D(2, 6) .* nds_3_1 - D(6, 6) .* nds_2_2 + D(6, 6) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(2, 1) .* nds_2_1 - D(2, 1) .* nds_3_1 - D(6, 1) .* nds_2_2 + D(6, 1) .* nds_3_2)) ./ (2 .* J);
el_K(:, 2, 8) = ((nds_1_1 - nds_2_1) .* (D(2, 2) .* nds_2_1 - D(2, 2) .* nds_3_1 - D(6, 2) .* nds_2_2 + D(6, 2) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(2, 6) .* nds_2_1 - D(2, 6) .* nds_3_1 - D(6, 6) .* nds_2_2 + D(6, 6) .* nds_3_2)) ./ (2 .* J);
el_K(:, 2, 9) = ((nds_1_1 - nds_2_1) .* (D(2, 4) .* nds_2_1 - D(2, 4) .* nds_3_1 - D(6, 4) .* nds_2_2 + D(6, 4) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(2, 5) .* nds_2_1 - D(2, 5) .* nds_3_1 - D(6, 5) .* nds_2_2 + D(6, 5) .* nds_3_2)) ./ (2 .* J);
el_K(:, 3, 1) = (((D(4, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(4, 1) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 1) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 3, 2) = (((D(4, 2) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 2) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(4, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 3, 3) = (((D(4, 4) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 4) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(4, 5) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 5) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 3, 4) = (((D(4, 1) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 1) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(4, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 3, 5) = (((D(4, 6) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(4, 2) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 2) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 3, 6) = (((D(4, 5) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 5) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(4, 4) .* (nds_2_1 - nds_3_1)) ./ J - (D(5, 4) .* (nds_2_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 3, 7) = ((nds_1_1 - nds_2_1) .* (D(4, 6) .* nds_2_1 - D(4, 6) .* nds_3_1 - D(5, 6) .* nds_2_2 + D(5, 6) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(4, 1) .* nds_2_1 - D(4, 1) .* nds_3_1 - D(5, 1) .* nds_2_2 + D(5, 1) .* nds_3_2)) ./ (2 .* J);
el_K(:, 3, 8) = ((nds_1_1 - nds_2_1) .* (D(4, 2) .* nds_2_1 - D(4, 2) .* nds_3_1 - D(5, 2) .* nds_2_2 + D(5, 2) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(4, 6) .* nds_2_1 - D(4, 6) .* nds_3_1 - D(5, 6) .* nds_2_2 + D(5, 6) .* nds_3_2)) ./ (2 .* J);
el_K(:, 3, 9) = ((nds_1_1 - nds_2_1) .* (D(4, 4) .* nds_2_1 - D(4, 4) .* nds_3_1 - D(5, 4) .* nds_2_2 + D(5, 4) .* nds_3_2)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(4, 5) .* nds_2_1 - D(4, 5) .* nds_3_1 - D(5, 5) .* nds_2_2 + D(5, 5) .* nds_3_2)) ./ (2 .* J);
el_K(:, 4, 1) = (((D(1, 6) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(1, 1) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 1) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 4, 2) = (((D(1, 2) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 2) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(1, 6) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 4, 3) = (((D(1, 4) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 4) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2 - (((D(1, 5) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 5) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2;
el_K(:, 4, 4) = (((D(1, 1) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 1) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(1, 6) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 4, 5) = (((D(1, 6) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 6) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(1, 2) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 2) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 4, 6) = (((D(1, 5) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 5) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2 - (((D(1, 4) .* (nds_1_2 - nds_3_2)) ./ J - (D(6, 4) .* (nds_1_1 - nds_3_1)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2;
el_K(:, 4, 7) = ((nds_1_1 - nds_2_1) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_3_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_3_1)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(1, 1) .* nds_1_2 - D(1, 1) .* nds_3_2 - D(6, 1) .* nds_1_1 + D(6, 1) .* nds_3_1)) ./ (2 .* J);
el_K(:, 4, 8) = ((nds_1_1 - nds_2_1) .* (D(1, 2) .* nds_1_2 - D(1, 2) .* nds_3_2 - D(6, 2) .* nds_1_1 + D(6, 2) .* nds_3_1)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_3_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_3_1)) ./ (2 .* J);
el_K(:, 4, 9) = ((nds_1_1 - nds_2_1) .* (D(1, 4) .* nds_1_2 - D(1, 4) .* nds_3_2 - D(6, 4) .* nds_1_1 + D(6, 4) .* nds_3_1)) ./ (2 .* J) - ((nds_1_2 - nds_2_2) .* (D(1, 5) .* nds_1_2 - D(1, 5) .* nds_3_2 - D(6, 5) .* nds_1_1 + D(6, 5) .* nds_3_1)) ./ (2 .* J);
el_K(:, 5, 1) = (((D(2, 1) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 1) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(2, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 5, 2) = (((D(2, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(2, 2) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 2) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 5, 3) = (((D(2, 5) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 5) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(2, 4) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 4) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 5, 4) = (((D(2, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(2, 1) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 1) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 5, 5) = (((D(2, 2) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 2) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(2, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 5, 6) = (((D(2, 4) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 4) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(2, 5) .* (nds_1_1 - nds_3_1)) ./ J - (D(6, 5) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 5, 7) = ((nds_1_2 - nds_2_2) .* (D(2, 1) .* nds_1_1 - D(2, 1) .* nds_3_1 - D(6, 1) .* nds_1_2 + D(6, 1) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_3_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_3_2)) ./ (2 .* J);
el_K(:, 5, 8) = ((nds_1_2 - nds_2_2) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_3_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(2, 2) .* nds_1_1 - D(2, 2) .* nds_3_1 - D(6, 2) .* nds_1_2 + D(6, 2) .* nds_3_2)) ./ (2 .* J);
el_K(:, 5, 9) = ((nds_1_2 - nds_2_2) .* (D(2, 5) .* nds_1_1 - D(2, 5) .* nds_3_1 - D(6, 5) .* nds_1_2 + D(6, 5) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(2, 4) .* nds_1_1 - D(2, 4) .* nds_3_1 - D(6, 4) .* nds_1_2 + D(6, 4) .* nds_3_2)) ./ (2 .* J);
el_K(:, 6, 1) = (((D(4, 1) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 1) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(4, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 6, 2) = (((D(4, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(4, 2) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 2) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 6, 3) = (((D(4, 5) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 5) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_2 - nds_3_2)) ./ 2 - (((D(4, 4) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 4) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_2_1 - nds_3_1)) ./ 2;
el_K(:, 6, 4) = (((D(4, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(4, 1) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 1) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 6, 5) = (((D(4, 2) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 2) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(4, 6) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 6) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 6, 6) = (((D(4, 4) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 4) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_1 - nds_3_1)) ./ 2 - (((D(4, 5) .* (nds_1_1 - nds_3_1)) ./ J - (D(5, 5) .* (nds_1_2 - nds_3_2)) ./ J) .* (nds_1_2 - nds_3_2)) ./ 2;
el_K(:, 6, 7) = ((nds_1_2 - nds_2_2) .* (D(4, 1) .* nds_1_1 - D(5, 1) .* nds_1_2 - D(4, 1) .* nds_3_1 + D(5, 1) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(4, 6) .* nds_1_1 - D(5, 6) .* nds_1_2 - D(4, 6) .* nds_3_1 + D(5, 6) .* nds_3_2)) ./ (2 .* J);
el_K(:, 6, 8) = ((nds_1_2 - nds_2_2) .* (D(4, 6) .* nds_1_1 - D(5, 6) .* nds_1_2 - D(4, 6) .* nds_3_1 + D(5, 6) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(4, 2) .* nds_1_1 - D(5, 2) .* nds_1_2 - D(4, 2) .* nds_3_1 + D(5, 2) .* nds_3_2)) ./ (2 .* J);
el_K(:, 6, 9) = ((nds_1_2 - nds_2_2) .* (D(4, 5) .* nds_1_1 - D(5, 5) .* nds_1_2 - D(4, 5) .* nds_3_1 + D(5, 5) .* nds_3_2)) ./ (2 .* J) - ((nds_1_1 - nds_2_1) .* (D(4, 4) .* nds_1_1 - D(5, 4) .* nds_1_2 - D(4, 4) .* nds_3_1 + D(5, 4) .* nds_3_2)) ./ (2 .* J);
el_K(:, 7, 1) = ((nds_2_2 - nds_3_2) .* (D(1, 1) .* nds_1_2 - D(1, 1) .* nds_2_2 - D(6, 1) .* nds_1_1 + D(6, 1) .* nds_2_1)) ./ (2 .* J) - ((nds_2_1 - nds_3_1) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_2_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 2) = ((nds_2_2 - nds_3_2) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_2_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_2_1)) ./ (2 .* J) - ((nds_2_1 - nds_3_1) .* (D(1, 2) .* nds_1_2 - D(1, 2) .* nds_2_2 - D(6, 2) .* nds_1_1 + D(6, 2) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 3) = ((nds_2_2 - nds_3_2) .* (D(1, 5) .* nds_1_2 - D(1, 5) .* nds_2_2 - D(6, 5) .* nds_1_1 + D(6, 5) .* nds_2_1)) ./ (2 .* J) - ((nds_2_1 - nds_3_1) .* (D(1, 4) .* nds_1_2 - D(1, 4) .* nds_2_2 - D(6, 4) .* nds_1_1 + D(6, 4) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 4) = ((nds_1_1 - nds_3_1) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_2_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_2_1)) ./ (2 .* J) - ((nds_1_2 - nds_3_2) .* (D(1, 1) .* nds_1_2 - D(1, 1) .* nds_2_2 - D(6, 1) .* nds_1_1 + D(6, 1) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 5) = ((nds_1_1 - nds_3_1) .* (D(1, 2) .* nds_1_2 - D(1, 2) .* nds_2_2 - D(6, 2) .* nds_1_1 + D(6, 2) .* nds_2_1)) ./ (2 .* J) - ((nds_1_2 - nds_3_2) .* (D(1, 6) .* nds_1_2 - D(1, 6) .* nds_2_2 - D(6, 6) .* nds_1_1 + D(6, 6) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 6) = ((nds_1_1 - nds_3_1) .* (D(1, 4) .* nds_1_2 - D(1, 4) .* nds_2_2 - D(6, 4) .* nds_1_1 + D(6, 4) .* nds_2_1)) ./ (2 .* J) - ((nds_1_2 - nds_3_2) .* (D(1, 5) .* nds_1_2 - D(1, 5) .* nds_2_2 - D(6, 5) .* nds_1_1 + D(6, 5) .* nds_2_1)) ./ (2 .* J);
el_K(:, 7, 7) = (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(1, 1) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 1) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2 - (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(1, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2;
el_K(:, 7, 8) = (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(1, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2 - (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(1, 2) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 2) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2;
el_K(:, 7, 9) = (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(1, 5) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 5) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2 - (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(1, 4) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) - D(6, 4) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J))) ./ 2;
el_K(:, 8, 1) = ((nds_2_1 - nds_3_1) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_2_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(2, 1) .* nds_1_1 - D(2, 1) .* nds_2_1 - D(6, 1) .* nds_1_2 + D(6, 1) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 2) = ((nds_2_1 - nds_3_1) .* (D(2, 2) .* nds_1_1 - D(2, 2) .* nds_2_1 - D(6, 2) .* nds_1_2 + D(6, 2) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_2_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 3) = ((nds_2_1 - nds_3_1) .* (D(2, 4) .* nds_1_1 - D(2, 4) .* nds_2_1 - D(6, 4) .* nds_1_2 + D(6, 4) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(2, 5) .* nds_1_1 - D(2, 5) .* nds_2_1 - D(6, 5) .* nds_1_2 + D(6, 5) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 4) = ((nds_1_2 - nds_3_2) .* (D(2, 1) .* nds_1_1 - D(2, 1) .* nds_2_1 - D(6, 1) .* nds_1_2 + D(6, 1) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_2_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 5) = ((nds_1_2 - nds_3_2) .* (D(2, 6) .* nds_1_1 - D(2, 6) .* nds_2_1 - D(6, 6) .* nds_1_2 + D(6, 6) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(2, 2) .* nds_1_1 - D(2, 2) .* nds_2_1 - D(6, 2) .* nds_1_2 + D(6, 2) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 6) = ((nds_1_2 - nds_3_2) .* (D(2, 5) .* nds_1_1 - D(2, 5) .* nds_2_1 - D(6, 5) .* nds_1_2 + D(6, 5) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(2, 4) .* nds_1_1 - D(2, 4) .* nds_2_1 - D(6, 4) .* nds_1_2 + D(6, 4) .* nds_2_2)) ./ (2 .* J);
el_K(:, 8, 7) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(2, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(2, 1) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 1) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;
el_K(:, 8, 8) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(2, 2) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 2) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(2, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;
el_K(:, 8, 9) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(2, 4) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 4) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(2, 5) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(6, 5) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;
el_K(:, 9, 1) = ((nds_2_1 - nds_3_1) .* (D(4, 6) .* nds_1_1 - D(4, 6) .* nds_2_1 - D(5, 6) .* nds_1_2 + D(5, 6) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(4, 1) .* nds_1_1 - D(4, 1) .* nds_2_1 - D(5, 1) .* nds_1_2 + D(5, 1) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 2) = ((nds_2_1 - nds_3_1) .* (D(4, 2) .* nds_1_1 - D(4, 2) .* nds_2_1 - D(5, 2) .* nds_1_2 + D(5, 2) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(4, 6) .* nds_1_1 - D(4, 6) .* nds_2_1 - D(5, 6) .* nds_1_2 + D(5, 6) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 3) = ((nds_2_1 - nds_3_1) .* (D(4, 4) .* nds_1_1 - D(4, 4) .* nds_2_1 - D(5, 4) .* nds_1_2 + D(5, 4) .* nds_2_2)) ./ (2 .* J) - ((nds_2_2 - nds_3_2) .* (D(4, 5) .* nds_1_1 - D(4, 5) .* nds_2_1 - D(5, 5) .* nds_1_2 + D(5, 5) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 4) = ((nds_1_2 - nds_3_2) .* (D(4, 1) .* nds_1_1 - D(4, 1) .* nds_2_1 - D(5, 1) .* nds_1_2 + D(5, 1) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(4, 6) .* nds_1_1 - D(4, 6) .* nds_2_1 - D(5, 6) .* nds_1_2 + D(5, 6) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 5) = ((nds_1_2 - nds_3_2) .* (D(4, 6) .* nds_1_1 - D(4, 6) .* nds_2_1 - D(5, 6) .* nds_1_2 + D(5, 6) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(4, 2) .* nds_1_1 - D(4, 2) .* nds_2_1 - D(5, 2) .* nds_1_2 + D(5, 2) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 6) = ((nds_1_2 - nds_3_2) .* (D(4, 5) .* nds_1_1 - D(4, 5) .* nds_2_1 - D(5, 5) .* nds_1_2 + D(5, 5) .* nds_2_2)) ./ (2 .* J) - ((nds_1_1 - nds_3_1) .* (D(4, 4) .* nds_1_1 - D(4, 4) .* nds_2_1 - D(5, 4) .* nds_1_2 + D(5, 4) .* nds_2_2)) ./ (2 .* J);
el_K(:, 9, 7) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(4, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(4, 1) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 1) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;
el_K(:, 9, 8) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(4, 2) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 2) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(4, 6) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 6) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;
el_K(:, 9, 9) = (J .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) .* (D(4, 4) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 4) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2 - (J .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J) .* (D(4, 5) .* ((nds_1_1 - nds_3_1) ./ J - (nds_2_1 - nds_3_1) ./ J) - D(5, 5) .* ((nds_1_2 - nds_3_2) ./ J - (nds_2_2 - nds_3_2) ./ J))) ./ 2;

%Mass matrix
el_M = zeros(size(els, 1), 9, 9);
el_M(:, 1, 1) = (J .* rho) ./ 6;
el_M(:, 2, 2) = (J .* rho) ./ 6;
el_M(:, 3, 3) = (J .* rho) ./ 6;
el_M(:, 4, 4) = (J .* rho) ./ 6;
el_M(:, 5, 5) = (J .* rho) ./ 6;
el_M(:, 6, 6) = (J .* rho) ./ 6;
el_M(:, 7, 7) = (J .* rho) ./ 6;
el_M(:, 8, 8) = (J .* rho) ./ 6;
el_M(:, 9, 9) = (J .* rho) ./ 6;


%Damping matrix (Rayleigh)
alpha = rayleigh_coefs(1);
beta = rayleigh_coefs(2);
el_C = zeros(size(els, 1), 9, 9);
el_C(:, 1, 1) = alpha*el_M(:, 1, 1) + beta*el_K(:, 1, 1);
el_C(:, 2, 2) = alpha*el_M(:, 2, 2) + beta*el_K(:, 2, 2);
el_C(:, 3, 3) = alpha*el_M(:, 3, 3) + beta*el_K(:, 3, 3);
el_C(:, 4, 4) = alpha*el_M(:, 4, 4) + beta*el_K(:, 4, 4);
el_C(:, 5, 5) = alpha*el_M(:, 5, 5) + beta*el_K(:, 5, 5);
el_C(:, 6, 6) = alpha*el_M(:, 6, 6) + beta*el_K(:, 6, 6);
el_C(:, 7, 7) = alpha*el_M(:, 7, 7) + beta*el_K(:, 7, 7);
el_C(:, 8, 8) = alpha*el_M(:, 8, 8) + beta*el_K(:, 8, 8);
el_C(:, 9, 9) = alpha*el_M(:, 9, 9) + beta*el_K(:, 9, 9);

%CRemove unwanted DOFs from element matrices
[loc_nd, loc_df, el_K, el_C, el_M] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, el_K, el_C, el_M);

end
