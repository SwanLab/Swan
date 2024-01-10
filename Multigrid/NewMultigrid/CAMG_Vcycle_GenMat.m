function [A_list, P_list, max_level] = CAMG_Vcycle_GenMat(A, direct_n)
% Generate the coefficient matrices for Classic AMG V-cycle
% A         : [in]  Original coefficient matrix
% direct_n  : [in]  Threshold for the smallest size of A in V-cycle
% A_list    : [out] The array of coefficient matrices on each level
% P_list    : [out] The array of interpolation (prolongation) operators on each level
% max_level : [out] Maximum level for V-cycle, initial is 1
	n = size(A, 1);
	A_list = {};
	P_list = {};
	level  = 1;
	A_list(level) = {A};
	
	while (n > direct_n)
		[CA, P] = CAMG_GenerateCoarseProblem(A);
		P_list(level) = {P};
		A_list(level + 1) = {CA};
		level = level + 1;
		A = CA;
		n = size(A, 1);
		
		nnz_ratio = nnz(A) / (size(A, 1) * size(A, 2));
		fprintf('Level %d, matrix size = %d,\t non-zero percentage = %2.4f \n', level, size(A, 1), nnz_ratio * 100.0);
	end
	max_level = level;
end