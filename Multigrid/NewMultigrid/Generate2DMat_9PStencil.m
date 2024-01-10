function A = Generate2DMat_9PStencil(ny, nx)
% Generate a ny rows * nx columns 2D grid and the coefficient matrix using 
% 9-point stencil. Index for (iy, ix) is (iy - 1) * nx + ix
	n = ny * nx;
	nnz = 0;
	rows = zeros(9 * n, 1); 
	cols = zeros(9 * n, 1); 
	vals = zeros(9 * n, 1); 
	
	for iy = 1 : ny
		for ix = 1 : nx
			curr_pos = (iy - 1) * nx + ix;
			for kx = -1 : 1
				for ky = -1 : 1
					new_pos = curr_pos + ky * nx + kx;
					new_ix = ix + kx;
					new_iy = iy + ky;
					if ((new_ix >= 1) && (new_ix <= nx) && (new_iy >= 1) && (new_iy <= ny))
						nnz = nnz + 1;
						rows(nnz) = curr_pos;
						cols(nnz) = new_pos;
						vals(nnz) = -1;
						if ((kx == 0) && (ky == 0))
							vals(nnz) = 8;
						end
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