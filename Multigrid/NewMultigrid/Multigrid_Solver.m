function [x, vcycle_cnt, rel_res] = Multigrid_Solver(A, b, smoother, pre_steps, pos_steps, rn_tol, max_vcycle)
% Multigrid solver for A * x = b using Classic AMG
% A          : [in]  The inital coefficient matrix
% b          : [in]  The right hand side
% smoother   : [in]  Function handle for a iterative method as a smoother
% pre_steps  : [in]  Number of iterations in the pre-smoothing
% pos_steps  : [in]  Number of iterations in the post-smoothing
% rn_tol     : [in]  The tolerance of the relative residual norm
% max_vcycle : [in]  Maximum number of V-clcyes to perform
% x          : [out] The solution after vcycle_cnt V-clcye
% vcycle_cnt : [out] Numbers of V-clcyes performed
% rel_res    : [out] Relative residual norm after each V-clcye

	if (nargin < 3) smoother   = @GS_Iter; end
	if (nargin < 4) pre_steps  = 1;	       end
	if (nargin < 5) pos_steps  = 1;	       end
	if (nargin < 6) rn_tol     = 1e-10;    end
	if (nargin < 7) max_vcycle = 100;      end
	
	n  = size(A, 1);
	x  = zeros(n, 1);
	rn = norm(b);
	vcycle_cnt = 0;
	rel_res(1) = rn;
	rn_stop    = rn * rn_tol;
	direct_n   = 16;
	PR_coef    = 1;
	
	% Generate coefficient matrices and interpolation operators of each level at once
	tic;
	[A_list, P_list, max_level] = CAMG_Vcycle_GenMat(A, direct_n);
	gm_t = toc;
	
	% Repeat V-cycle until converge
	tic;
	while ((rn > rn_stop) && (vcycle_cnt <= max_vcycle))
		x  = Multigrid_Vcycle(1, A_list, P_list, b, x, direct_n, PR_coef, smoother, pre_steps, pos_steps);
		r  = b - A * x;
		rn = norm(r, 2);
		vcycle_cnt = vcycle_cnt + 1;
		rel_res(vcycle_cnt + 1) = rn / rel_res(1);
	end
	rel_res = rel_res(2 : end);
	vcyc_t = toc;
	
	fprintf('Matrices generating wall-time = %f (s)\n', gm_t);
	fprintf('V-cycle solver      wall-time = %f (s)\n', vcyc_t);
	fprintf('Performed V-cycles = %d\n', vcycle_cnt);
	fprintf('||b - A * x||_2    = %e\n', rn);
end