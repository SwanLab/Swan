function [CA, P] = CAMG_GenerateCoarseProblem(A)
% Using Classic AMG to generate the coarse grid matrix 
% and interlopation operator according to the input matrix
% A  : [in]  Fine grid problem matrix
% CA : [out] Coarse grid problem coefficient matrix
% P  : [out] Interpolation operator, CA = P^T * A * P
	
	is_cg_point = CAMG_GenerateCoarseGrid(A);
	cg_point_id = cumsum(is_cg_point);
	CA_n = sum(is_cg_point);
	A_n  = size(A, 1);
	sd_A = CAMG_GetStrongDependencyMat(A, 0.5);
	
	% Construct interpolation operator P
	P_rows = []; P_cols = []; P_vals = []; P_nnz = 0;
	n = size(A, 1);
	for curr_pid = 1 : n
		% P is CA_n cols, A_n rows, each row is a interplation relationship 
		% from coarse grid to fine grid
		[ip_id, ip_coef] = CAMG_InterpolateForPoint(curr_pid, A, sd_A, is_cg_point);
		for i = 1 : size(ip_id, 1)
			P_nnz = P_nnz + 1;
			P_rows(P_nnz) = curr_pid;
			P_cols(P_nnz) = cg_point_id(ip_id(i));
			P_vals(P_nnz) = ip_coef(i);
		end
	end
	P = sparse(P_rows, P_cols, P_vals, A_n, CA_n, P_nnz);
	
	% Construct coarse grid matrix
	CA = P' * A * P;
end