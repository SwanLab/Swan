classdef FEM < handle
    %FEM Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        problemID
        geometry
        mesh
        dof
        element
        Fext
        LHS
        variables
    end
    
    %% Restricted properties definition ===================================
    properties (GetAccess = ?Postprocess, SetAccess = private)
    end
    
    %% Private properties definition ======================================
    properties (Access = protected)
        solver
    end
    
    %% Public methods definition ==========================================
    methods (Static, Access = public)
        function obj = create(problemID)
            mesh = Mesh(problemID); % Mesh defined twice, but almost free
            switch mesh.ptype
                case 'ELASTIC'
                    obj = Elastic_Problem(problemID);
                case 'THERMAL'
                    obj = Thermal_Problem(problemID);
                case 'DIFF-REACT'
                    obj = DiffReact_Problem(problemID);
                case 'HYPERELASTIC'
                    obj = Hyperelastic_Problem(problemID);
                case 'Stokes'
                    obj = Stokes_Problem(problemID);
            end
        end
    end
    
    methods
        function print(obj)
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,results);
        end
        
        % !! Remove residual & change to more straight-forward sol?? !!
        function sol = solve_steady_problem(obj,free_dof)
            total_free_dof = sum(free_dof);
            dr = obj.element.computedr;
            x0 = zeros(total_free_dof,1);
            r = obj.element.computeResidual(x0,dr);
            sol = obj.solver.solve(dr,-r);
        end
        
        function setDof(obj,dof)
            obj.dof = dof;
        end
        
        function setMatProps(obj,props)
            obj.element.material = obj.material.setProps(props);
        end
    end
    
    methods (Abstract)
        preProcess(obj)
        computeVariables(obj)
        postProcess(obj)
    end 
end
