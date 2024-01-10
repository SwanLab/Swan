function [ip_id, ip_coef] = CAMG_InterpolateForPoint(curr_pid, A, sd_A, is_cg)
% Get the interpolation relationship of a given point
% curr_pid : [in]  ID of grid point need to be interpolated
% A        : [in]  Original coefficient matrix
% sd_A     : [in]  Strong dependency matrix, from CAMG_GetStrongDependencyMat()
% is_cg    : [in]  If a grid point is a coarse grid point, from CAMG_GenerateCoarseGrid()
% ip_id    : [out] IDs of coarse grid points used to interpolate 
% ip_coef  : [out] Weights of coarse grid points used to interpolate 
	
	% If current point is coarse grid point, need not to interpolate
	if (is_cg(curr_pid) == 1)
		ip_id   = curr_pid;
		ip_coef = 1;
		return;
	end
	
	% Find all coarse grid points for interpolation
	neighbors = find(A(curr_pid, :) ~= 0);
	neighbors = setdiff(neighbors, curr_pid); 
	neighbors_cnt = size(neighbors, 2);
	ip_cnt  = 0;
	ip_id   = zeros(neighbors_cnt, 1);
	ip_coef = zeros(neighbors_cnt, 1);
	for i = 1 : neighbors_cnt
		curr_neighbor = neighbors(i);
		if (is_cg(curr_neighbor) == 1)
			ip_cnt = ip_cnt + 1;
			ip_id(ip_cnt)   = curr_neighbor;
			ip_coef(ip_cnt) = A(curr_pid, curr_neighbor);
		end
	end
	ip_id = ip_id(1 : ip_cnt);
	ip_coef = ip_coef(1 : ip_cnt);
	
	% Redistribute the edge weight of all fine grid points with weak dependency
	a_ii_new = A(curr_pid, curr_pid);
	for i = 1 : neighbors_cnt
		curr_neighbor = neighbors(i);
		if ((is_cg(curr_neighbor) == 0) && (sd_A(curr_pid, curr_neighbor) == -1))
			a_ii_new = a_ii_new + A(curr_pid, curr_neighbor);
		end
	end
	
	% Redistribute the edge weight of all fine grid points with strong dependency
	for i = 1 : neighbors_cnt
		curr_neighbor = neighbors(i);
		if ((is_cg(curr_neighbor) == 0) && (sd_A(curr_pid, curr_neighbor) == 1))
			cg_neighbor_list = [];
			cg_neighbor_weight_sum = 0;
			
			% Find all coarse grid points connected to curr_neighbor
			% and sum the edge weight of these connections
			for j = 1 : neighbors_cnt
				curr_check_neighbor = neighbors(j);
				if ((is_cg(curr_check_neighbor) == 1) && (A(curr_neighbor, curr_check_neighbor) ~= 0))
					cg_neighbor_weight_sum = cg_neighbor_weight_sum + A(curr_neighbor, curr_check_neighbor);
					cg_neighbor_list = [cg_neighbor_list j];
				end
			end
			
			% Redistribute the edge weight A(curr_pid, curr_neighbor) to each 
			% neighboring coarse grid point
			for k = 1 : size(cg_neighbor_list, 2)
				j = cg_neighbor_list(k);
				curr_check_neighbor = neighbors(j);
				contrib = A(curr_neighbor, curr_check_neighbor) / cg_neighbor_weight_sum * A(curr_pid, curr_neighbor);
				
				ip_id_pos = find(ip_id == curr_check_neighbor);
				ip_coef(ip_id_pos) = ip_coef(ip_id_pos) + contrib;
			end
		end
	end
	
	% Normalize the coefficients
	ip_coef = ip_coef ./ -a_ii_new;
end