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
            obj.mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Stokes(problemID,obj.geometry);
        end
        
        function preProcess(obj)
            cParams.nelem = obj.geometry(1).interpolation.nelem;
            obj.material = Material_Stokes(cParams);
            obj.element = Element_Stokes(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end
            transient = false;  % !! This should not be defined in here !!
            tol = 1e-6;         % !! This should not be defined in here !! 
            if transient
                dt = 0.01;      % !! This should not be defined in here !!
                final_time = 1; % !! This should not be defined in here !!
                x = obj.solve_transient_nonlinear_problem(free_dof,tol,dt,final_time);
            else
                x = obj.solve_steady_nonlinear_problem(free_dof,tol);
            end
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
        
        function createGeometry(obj,mesh)
                obj.geometry = Geometry(mesh,'QUADRATIC');
                obj.geometry(2) = Geometry(mesh,'LINEAR');
                obj.geometry(1).nfields = 2;
        end
    end
end

