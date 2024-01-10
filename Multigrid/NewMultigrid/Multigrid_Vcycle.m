function x = Multigrid_Vcycle(level, A_list, P_list, b, x0, direct_n, PR_coef, smoother, pre_steps, pos_steps)
% Multigrid V-clcye
% level     : [in]  The current level (initial is 1)
% A_list    : [in]  The array of coefficient matrices on each level
% P_list    : [in]  The array of interpolation (prolongation) operators on each level
% b         : [in]  The right hand side
% x0        : [in]  Initial guess
% direct_n  : [in]  Threshold for directly solving A_list(level) * x = b
% PR_coef   : [in]  The coefficient constant between restriction and prolongation
% smoother  : [in]  Function handle for a iterative method as a smoother
% pre_steps : [in]  Number of iterations in the pre-smoothing
% pos_steps : [in]  Number of iterations in the post-smoothing
% x         : [out] The new solution after a V-clcye

	% Load coefficient matrix
	A = cell2mat(A_list(level));
	
	% If the problem is small enough, solve it directly
	n = size(b, 1);
	if (n <= direct_n)
		x = A \ b;
		return;
	end

	% Pre-smoothing
	x = smoother(A, b, 1e-14, pre_steps, x0);
	
	% Load restriction operator and construct interpolation operator
	P = cell2mat(P_list(level));
	R = P' * PR_coef;
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	r   = b - A * x;
	r_C = R * r;
	
	% Solve coarse grid problem recursively
	x0  = zeros(coarse_n, 1);
	e_C = Multigrid_Vcycle(level + 1, A_list, P_list, r_C, x0, direct_n, PR_coef, smoother, pre_steps, pos_steps);
	
	% Transfer error to fine grid and correct
	x = x + P * e_C;
	
	% Post-smoothing
	x = smoother(A, b, 1e-14, pos_steps, x);
end