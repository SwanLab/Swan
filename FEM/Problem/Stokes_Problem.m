classdef Stokes_Problem < FEM
    %Stokes_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    %% Restricted properties definition ===================================
    properties %(GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Stokes_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Stokes(problemID,obj.geometry);
        end
        
        function preProcess(obj)
            obj.material = Material_Stokes(obj.geometry(1).interpolation.nelem);
            obj.element = Element_Stokes(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            transient = false;  % !! This should not be defined in here !!
            if transient
                tol = 1e-6;     % !! This should not be defined in here !! 
                dt = 0.01;      % !! This should not be defined in here !!
                final_time = 1; % !! This should not be defined in here !!
                x = obj.solve_transient_problem(free_dof,tol,dt,final_time);
            else
                x = obj.solve_steady_problem(free_dof);
            end
            obj.variables = obj.element.computeVars(x);
        end
        
        function sol = solve_transient_problem(obj,free_dof,tol,dt,final_time)
            total_free_dof = sum(free_dof);
            x_n(:,1) = zeros(total_free_dof,1);
            x0 = zeros(total_free_dof,1);
            
            dr = obj.element.computedr(dt);
            
            for istep = 2: final_time/dt
                u_previous_step = x_n(1:free_dof(1),istep-1);
                
                r = obj.element.computeResidual(x0,dr,u_previous_step);
                while dot(r,r) > tol
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
        
%         function print(obj)
%             postprocess = Postprocess_PhysicalProblem();
%             results.physicalVars = obj.variables;
%             postprocess.print(obj,obj.problemID,results);
%         end
        
        function postProcess(obj)
            % ToDo
            % Inspire in TopOpt
            
        end
        
        function createGeometry(obj,mesh)
                obj.geometry = Geometry(mesh,'QUADRATIC');
                obj.geometry(2) = Geometry(mesh,'LINEAR');
                obj.geometry(1).nfields = 2;
        end
    end
end

