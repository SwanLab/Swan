function test_Poisson2D(n)
	N = n * n;        % Number of inner grid points on (0,1)*(0,1)
	rng(n);
	fprintf('Using %d * %d square initial grid\n', n, n);
	A = Generate2DMat_5PStencil(n, n);
	b = rand(N, 1) - 0.5;
	[x, vc_cnt, rel_res] = Multigrid_Solver(A, b, @GS_Iter, 1, 1, 1e-10, 20);
end