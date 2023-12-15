function v_s = smoothEll(obj, A, v, g, k, method)    
	
    if strcmp(method, 'Jacobi')	
        method_func = @(v_prev) obj.JacobiIter(A, v_prev, g);
    elseif strcmp(method, 'Seidel')
        method_func = @(v_prev) obj.SeidelIter(A, v_prev, g);
    elseif strcmp(method, 'SOR')
		omega_star = 1.9;
        method_func = @(v_prev) obj.SORIter(A, v_prev, g, omega_star);
    else
        disp('smooth::ERROR: undefined method')
        v_s = [];
        return;
    end
    
    v_s = v;
    for kk = 1:k     
        v_s = method_func(v_s);
    end 
    
end
