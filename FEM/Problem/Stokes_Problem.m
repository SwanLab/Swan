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
    
    properties (Access = private)
        interp
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Stokes_Problem(fileName)
            obj.readProblemData(fileName);
            obj.createGeometry(obj.mesh);
            obj.createInterpolation();
            obj.dof = DOF_Stokes(fileName,obj.mesh,obj.geometry,obj.interp);
            cParams.nelem = obj.mesh.nelem;
            obj.material = Material_Stokes(cParams);
            obj.element  = Element_Stokes(obj.geometry,obj.mesh,obj.material,obj.dof,obj.problemData,obj.interp);
            obj.solver   = Solver.create;
        end
        
        function computeVariables(obj)
            nFields = numel(obj.interp);
            for ifield = 1:nFields
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
            s.mesh = mesh;
            obj.geometry    = Geometry.create(s);
            obj.geometry(2) = Geometry.create(s);
            %obj.geometry(1).nfields = 2;
        end
    end
    
    methods (Access = private)
        
        function createInterpolation(obj)
            interpU = 'QUADRATIC';
            interpP = 'LINEAR';
            obj.interp{1}=Interpolation.create(obj.mesh,interpU);
            obj.interp{2}=Interpolation.create(obj.mesh,interpP);
        end
    end
end

