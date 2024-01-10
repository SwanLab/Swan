function sd_mat = CAMG_GetStrongDependencyMat(A, theta)
% Mark all strong dependency in the graph that 
% $-a_{ij} >= theta * max_{k != i} {-a_{ik}}$
% A      : [in]  Edge list of a weighted directed graph 
% theta  : [in]  Strong dependency threshold
% sd_mat : [out] Edge list, marked strong connections: 1 for true, -1 for false
	n = size(A, 1);
	sd_mat = sparse(n, n);
	for row = 1 : n
		threshold  = theta * max(-A(row, :));
		if (threshold == 0)
			continue;
		end
		nnz_in_row = find( A(row, :) ~= 0);         
		nnz_in_row = setdiff(nnz_in_row, [row]);    % Ignore diagonal element
		sd_in_row  = find(-A(row, :) >= threshold);
		wd_in_row  = setdiff(nnz_in_row, sd_in_row);
		sd_mat(row, sd_in_row) =  1;
		sd_mat(row, wd_in_row) = -1;
	end
end