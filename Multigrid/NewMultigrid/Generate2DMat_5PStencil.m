function A = Generate2DMat_9PStencil(ny, nx)
% Generate a ny rows * nx columns 2D grid and the coefficient matrix using 
% 5-point stencil. Index for (iy, ix) is (iy - 1) * nx + ix
	n = ny * nx;
	nnz = 0;
	rows = zeros(9 * n, 1); 
	cols = zeros(9 * n, 1); 
	vals = zeros(9 * n, 1); 
	dx = [0, 0,  0, 1, -1];
	dy = [0, 1, -1, 0,  0];
	
	for iy = 1 : ny
		for ix = 1 : nx
			curr_pos = (iy - 1) * nx + ix;
			for k = 1 : 5
				new_ix = ix + dx(k);
				new_iy = iy + dy(k);
				new_pos = (new_iy - 1) * nx + new_ix;
				if ((new_ix >= 1) && (new_ix <= nx) && (new_iy >= 1) && (new_iy <= ny))
					nnz = nnz + 1;
					rows(nnz) = curr_pos;
					cols(nnz) = new_pos;
					vals(nnz) = -1;
					if (k == 1)
						vals(nnz) = 4;
					end
				end
			end
		end
	end
	
	rows = rows(1 : nnz);
	cols = cols(1 : nnz);
	vals = vals(1 : nnz);
	A = sparse(rows, cols, vals, n, n, nnz);
end