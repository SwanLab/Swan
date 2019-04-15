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
        function setupFromGiDFile(obj,problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
            obj.mesh.ptype = 'DIFF-REACT';
        end
        
        function setupFromMesh(obj,mesh)
            obj.mesh = mesh; % Mesh defined twice, but almost free
            obj.problemID = mesh.problemID;
        end
        
        function preProcess(obj)
            obj.createGeometry(obj.mesh);
            obj.setDOFs();
            % !! Material required?? !!
            obj.setElement;
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj,x)
            bc = obj.element.getBcApplier();
            x_red  = bc.full_vector_2_reduced_vector(x);
            LHS = obj.element.computeLHS;
            x_reg = obj.solver.solve(LHS,x_red);
            obj.variables.x = bc.reduced_vector_2_full_vector(x_reg);
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
    
    methods (Access = protected)
        function setElement(obj)
            obj.element = Element_DiffReact(obj.mesh,obj.geometry,obj.material,obj.dof);
        end
        
        function setDOFs(obj)
            obj.dof = DOF_DiffReact(obj.geometry);
        end
    end
end
