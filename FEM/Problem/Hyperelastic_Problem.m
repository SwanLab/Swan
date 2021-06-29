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
            obj.mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
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
            x = obj.solve_steady_nonlinear_problem(free_dof,tol);
            obj.variables = obj.element.computeVars(x);
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
        
        % !! THIS SHOULD BE DEFINED BY THE USER !!
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
    end
end

