classdef DiffReact_Problem < FEM
    %DiffReact_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    %% Restricted properties definition ===================================
    properties %(GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material % !! Preguntar si necessita material
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = DiffReact_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.mesh.ptype = 'DIFF-REACT';
            obj.dof = DOF_DiffReact(obj.geometry);
        end
        
        function preProcess(obj)
            % !! Material required?? !!
            obj.element = Element_DiffReact(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj,x)
            x_red  = obj.element.full_vector_2_reduced_vector(x);

            tol = 1e-6;
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
                obj.dof.neumann = obj.dof.free{ifield};
            end
            obj.dof.neumann_values = x_red;
            x_red = obj.solve_steady_problem(free_dof,tol);
            obj.variables.x = obj.element.reduced_vector_2_full_vector(x_red);
        end
        
        function postProcess(obj)
            % ToDo
            % Inspire in TopOpt
        end
        
        function obj = setEpsilon(obj,epsilon)
            obj.element.setEpsilon(epsilon);
        end
        
        % !! THIS SHOULD BE DEFINED BY THE USER !!
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
    end
end

