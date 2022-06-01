classdef NonLinear_Solver < handle

    properties (Access = public)
    end

    properties (Access = private)
        tol
        state
        solver
        element
        free_dof
    end

    methods (Access = public)

        function obj = NonLinear_Solver(cParams)
            obj.init(cParams);
            obj.createSolver();
        end

        function x = solve(obj, params)
            switch params.state
                case 'Steady'
                    x = obj.solveSteady();
                case 'Transient'
                    x = obj.solveTransient(params);
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.tol      = cParams.tol;
%             obj.state    = cParams.state;
            obj.element  = cParams.element;
            obj.free_dof = cParams.free_dof;
        end

        function createSolver(obj)
            s.type = 'DIRECT';
            obj.solver = Solver.create(s);
        end

        function sol = solveSteady(obj)
            total_free_dof = sum(obj.free_dof);
            dt = Inf;
            dr = obj.element.computedr(dt);
            x0 = zeros(total_free_dof,1);
            
            r = obj.element.computeResidual(x0,dr);
            x = obj.convergeSolution(dr, r, x0);
            sol = x;
        end
        
        function sol = solveTransient(obj,params)
            dt = params.dt;
            final_time = params.final_time;
            total_free_dof = sum(obj.free_dof);
            x_n(:,1) = zeros(total_free_dof,1);
            x0 = zeros(total_free_dof,1);
            
            dr = obj.element.computedr(dt);
            
            for istep = 2: final_time/dt
                u_previous_step = x_n(1:obj.free_dof(1),istep-1);
                
                r = obj.element.computeResidual(x0,dr,u_previous_step);
                while dot(r,r) > obj.tol
                    inc_x = obj.solver.solve(dr,-r);
                    x = x0 + inc_x;
                    % Compute r
                    r = obj.element.computeResidual(x,dr,u_previous_step);
                    x0 = x;
                end
                x_n(:,istep) = x;
            end
            sol = x_n;
        end

        function sol = convergeSolution(obj, dr, r, x0, u_previous_step)
            while dot(r,r) > obj.tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r
                r = obj.element.computeResidual(x,dr);
                x0 = x;
            end
            sol = x0;
        end

    end

end

