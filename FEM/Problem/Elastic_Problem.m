classdef Elastic_Problem < FEM
    %Elastic_Problem Summary of this class goes here
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
        function obj = Elastic_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Elastic(problemID,obj.geometry,obj.mesh);
            obj.material = Material.create(obj.geometry,obj.mesh);
        end
        
        function preProcess(obj)
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            Kred = obj.element.computeLHS;
            fext_red = obj.element.computeRHS;
            u = obj.solver.solve(Kred,fext_red);
            obj.variables = obj.element.computeVars(u);
        end
        
%         function print(obj)
%             postprocess = Postprocess_PhysicalProblem;
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

