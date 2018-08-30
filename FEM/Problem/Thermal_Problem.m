classdef Thermal_Problem < FEM
    %Thermal_Problem Summary of this class goes here
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
        function obj = Thermal_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF(problemID,obj.geometry);
        end
        
        function preProcess(obj)
            obj.element = Element_Thermal(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            x = obj.solve_steady_problem(free_dof);
            obj.variables = obj.element.computeVars(x);
        end
        
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

