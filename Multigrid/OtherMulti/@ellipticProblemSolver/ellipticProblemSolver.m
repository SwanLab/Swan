% Class 'ellipticProblemSolver' implements a solver of elliptic problem,
% includes iteration methods: Jacobi, Seidel, SOR and multigrid
classdef ellipticProblemSolver < handle
    properties (Access = public)
        % grid
        grid = []
        
        % finite-difference operator
        A_op = [];
        
        % right side of equation
        f_func = [];               
        % border conditions
        border = [];        
        % exact solution (if it's known)
        u_ex = [];       
        
        % spectral radius
        rho = [];        
        % error
        error = [];   

        % max number of iterations
        k_max = [];
        % accuracy
        eps = []; 

        % params of MG method
        nu_1 = [];
        nu_2 = [];
        level = [];
        smooth_m = []; 
    end
    
    methods (Access = public)
        function obj = ellipticProblemSolver(c, e_l, N,...
                a_coeff, b_coeff, q_coeff,...
                f, border, u_ex)
            obj.grid = gridClass(c, e_l, N);
            obj.A_op = finDiffOpClass(a_coeff, b_coeff, q_coeff);
            obj.f_func = f;
            obj.border = border;     
            obj.u_ex = u_ex;
        end       
        
        % PARAM IN:
        %   * grid_m - object of gridClass
        % PARAM_OUT:
        %   * border_vec - a vector with values of solution in border
        function border_vec = fillBorder(obj, grid_m)
            N = grid_m.N;
            border_vec = zeros(N + 2);
            
            % left
            [x, y] = grid_m.leftEdge();
            border_vec(:, 1)   = obj.border(x, y);
            % right
            [x, y] = grid_m.rightEdge();
            border_vec(:, end) = obj.border(x, y);
            % bottom
            [x, y] = grid_m.bottomEdge();
            border_vec(1,   :) = obj.border(x, y);
            % top
            [x, y] = grid_m.topEdge();
            border_vec(end, :) = obj.border(x, y);               
        end
        
        % PARAM_OUT:
        %   * g_vec - right part of a finite difference equation
        function g_vec = getRightPart(obj)
            N = obj.grid.N;
            g_vec = obj.fillBorder(obj.grid);
            
            % indices of inner nodes
            i_in = 2:N+1;
            j_in = 2:N+1;            
            
            g_vec(j_in, i_in) = obj.f_func(obj.grid.x(i_in),...
                                           obj.grid.y(j_in));
        end        
        
        % PARAM IN:
        %   * nu_1, nu_2 - number of pre- and post-smoothing iterations,
        %   * level_max  - number of grids,
        %   * smooth_m   - smoothing method
        % PARAM_OUT:
        %   * v - vector of numerical solution, multigrid method
        function v = solveProblemByMG(obj, nu_1, nu_2, level_max, smooth_m)            
            A_mat = obj.A_op.getMatrix(obj.grid);
            g = obj.getRightPart();
            v_prev = obj.fillBorder(obj.grid);
            
            absolute_error_flag = ~isempty(obj.u_ex);
            if absolute_error_flag
                u0 = obj.u_ex(obj.grid.x, obj.grid.y);
            end
            
            obj.error = ones(obj.k_max, 1);
            obj.rho = zeros(obj.k_max, 1);
            k = 1;
            while k <= obj.k_max && obj.eps < obj.error(max(k-1 ,1))
                v = MultiGrid(obj, A_mat, v_prev, g, 1, nu_1, nu_2, 1,...
                              level_max, smooth_m);
                if k > 1
                    rho_numerator = norm(v - v_prev);
                    obj.rho(k) = rho_numerator/rho_denominator;
                end

                if absolute_error_flag
                    obj.error(k) = norm(u0 - v);
                else
                    obj.error(k) = norm(v_prev - v);
                end
                rho_denominator = norm(v - v_prev);
                v_prev = v;
                
                sprintf('k = %d, err = %.2e\n', k, obj.error(k))
                k = k + 1;
            end   
            obj.error = obj.error(1:k-1);
            obj.rho = obj.rho(1:k-1);
        end
        
        % PARAM IN:
        %   * method - iteration method,
        %   * omega  - parameter of SOR method
        % PARAM_OUT:
        %   * v - vector of numerical solution
        function v = solveProblemByIter(obj, method, omega)
            absolute_error_flag = ~isempty(obj.u_ex);
            if absolute_error_flag
                u0 = obj.u_ex(obj.grid.x, obj.grid.y);
            end
            
            A_mat = obj.A_op.getMatrix(obj.grid);
            v_prev = obj.fillBorder(obj.grid);
            g = obj.getRightPart();
            
            if strcmp(method, 'Jacobi')	
                method_func = @(v_prev) obj.JacobiIter(A_mat, v_prev, g);
            elseif strcmp(method, 'Seidel')
                method_func = @(v_prev) obj.SeidelIter(A_mat, v_prev, g);
            elseif strcmp(method, 'SOR')
                method_func = @(v_prev) obj.SORIter(A_mat, v_prev, g, omega);
            else
                disp('solveProblemByIter::ERROR: undefined iteration method')
                v = [];
                return
            end
                                    
            obj.error = ones(obj.k_max, 1);	
            obj.rho = zeros(obj.k_max, 1);
            k = 1;
            while k <= obj.k_max && obj.eps < obj.error(max(k-1 ,1))
                v = method_func(v_prev);
                if k > 1
                    rho_numerator = norm(v - v_prev);
                    obj.rho(k) = rho_numerator/rho_denominator;
                end
                if absolute_error_flag
                    obj.error(k) = norm(u0 - v);                
                else
                    obj.error(k) = norm(v_prev - v);                
                end
                rho_denominator = norm(v - v_prev);
                v_prev = v;
                if ~rem(k, 100)
                    sprintf('k = %d, err = %.2e\n', k, obj.error(k))
                end
                k = k + 1;
            end
            obj.error = obj.error(1:k-1);            
            obj.rho = obj.rho(1:k-1);
        end
        
        % PARAM IN:
        %   * eps   - accuracy of error,
        %   * k_max - max number of iterations
        function setParams(obj, eps, k_max)
            obj.eps = eps;
            obj.k_max = k_max;
        end
        
        % PARAM IN:
        %   * A - matrix of problem,
        %   * v - approximate solution,
        %   * g - right part of equation,
        %   * k - number of smoothing iterations,
        %   * method - method of smoothing
        % PARAM_OUT:
        %   * v_s - smoothing solution of the problem Av = g
        v_s = smoothEll(obj, A, v, g, k, method)
    end
    
    methods (Static)        
        d_h = prolongation(d_2h)
        
        d_2h = restriction(d_h)   
        
        v = JacobiIter(A, v0, g)
        
        v = SeidelIter(A, v0, g)
        
        v = SORIter(A, v0, g, omega)
    end
    
    methods (Access = private)        
        vm_new = MultiGrid(obj, A_m, vm_old, g_m, gamma, nu_1, nu_2,...
                           level, level_max, smooth_m)                  
    end
          
end