function is_coarse_grid_point = CAMG_GenerateCoarseGrid(A)
% Generate a "coarse grid" in Classic AMG with "fine grid" A
% is_coarse_grid_point[i] = 1 if grid point i is a coarse grid point
%                    !!! NOTICE !!! 
% The "Second pass" is not included, since "Researchers later found 
% that the second pass leads to high computational costs, and they 
% have largely abandoned it in favor of other approaches."
	n = size(A, 1);
	edge_list = (A ~= 0);
	edge_list = edge_list - diag(diag(edge_list));  % Remove the diagonal 
	x_degree  = edge_list * ones(n, 1);
	cg_point  = -ones(n, 1);
	
	unvisited_points = sum(cg_point == -1);
	while unvisited_points > 0
		% Find next coarse grid point
		[~, new_cg_point] = max(x_degree); 
		cg_point(new_cg_point) = 1;
		x_degree(new_cg_point) = 0; 
		
		% Find new fine grid point neighbors
		neighbors = find(edge_list(new_cg_point, :) == 1);
		new_neighbors = [];
		for i = 1 : size(neighbors, 2)
			if (cg_point(neighbors(i)) == -1)
				new_neighbors = [new_neighbors neighbors(i)];
				% Mark new neighbors
				cg_point(neighbors(i)) = 0;  
				x_degree(neighbors(i)) = 0; 
			end
		end
		
		% Mark new neighbors and add their unvisited neighbors' weight
		for i = 1 : size(new_neighbors, 2)
			curr_neighbor = new_neighbors(i);
			
			nn = find(edge_list(curr_neighbor, :) == 1);
			nn2 = [];
			for i = 1 : size(nn, 2)
				if (cg_point(nn(i)) == -1)      % Only update those unvisited neighbors
					nn2 = [nn2 nn(i)];
				end
			end
			
			for i = 1 : size(nn2, 2)
				x_degree(nn2(i)) = x_degree(nn2(i)) + 1;
			end
		end
		
		unvisited_points = sum(cg_point == -1);
	end
	
	is_coarse_grid_point = cg_point;
end