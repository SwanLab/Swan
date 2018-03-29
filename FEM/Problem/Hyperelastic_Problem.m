classdef Hyperelastic_Problem < FEM
    %Hyperelastic_Problem Summary of this class goes here
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
        function obj = Hyperelastic_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Elastic(problemID,obj.geometry,obj.mesh);
            obj.material = Material.create(obj.geometry,obj.mesh);
        end
        
        function preProcess(obj)
            obj.element = Element_Hyperelastic(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            tol = 1e-6; % !! This should not be defined in here !!
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            x = obj.solve_steady_problem(free_dof,tol);
            obj.variables = obj.element.computeVars(x);
        end
        
        function print(obj)
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,results);
        end
        
        function postProcess(obj)
            % ToDo
            % Inspire in TopOpt
            
        end
        
        function sol = solve_steady_problem(obj,free_dof,tol)
            total_free_dof = sum(free_dof);
            dr = obj.element.computedr;
            x0 = zeros(total_free_dof,1);
            
            r = obj.element.computeResidual(x0,dr);
            x = x0;
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r
                r = obj.element.computeResidual(x,dr);
                x0 = x;
            end
            sol = x;
        end
        
        % !! THIS SHOULD BE DEFINED BY THE USER !!
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
    end
end

