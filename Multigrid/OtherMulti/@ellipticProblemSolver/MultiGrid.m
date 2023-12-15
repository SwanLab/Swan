% grids is numerated from 1 (first grid is obj.grid) 
% to MAX_GRID_NUM = log2(N+1) (grid without inner nodes - 
% just 4 border nodes) so
% param 'level_max' should be < MAX_GRID_NUM
function vm_new = MultiGrid(obj, A_m, vm_old, g_m, gamma, nu_1, nu_2,...
                            level, level_max, smooth_m)    
    
    if level < level_max
        
        vm_line = obj.smoothEll(A_m, vm_old, g_m, nu_1, smooth_m);
        
        dm = g_m - A_m.mul(vm_line);
        dm_1 = obj.restriction(dm);
        
        em_1 = zeros(size(dm_1));        
        grid_m_1 = obj.grid.getSubgrid(level+1);
        Am_1 = obj.A_op.getMatrix(grid_m_1);
        for i = 1:gamma
            em_1 = obj.MultiGrid(Am_1, em_1, dm_1, gamma, nu_1, nu_2,...
                                 level+1, level_max, smooth_m);
        end
        
        em_wave = obj.prolongation(em_1);        
        vm_wave = vm_line + em_wave;
        vm_new = obj.smoothEll(A_m, vm_wave, g_m, nu_2, smooth_m);
        
    else
        vm_new = A_m.solveSystem(g_m);
    end    
end