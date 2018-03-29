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
        
        function print(obj,file_name)
            postprocess = Postprocess_PhysicalProblem;
            if nargin > 1
                postprocess.print(obj,file_name,obj.variables);
            else
                postprocess.print(obj,obj.problemID,obj.variables);
            end
        end
        
        function syncPostProcess(obj,evtobj)
            addlistener(evtobj,'res_file','PostSet',@obj.print_slave);
        end
        
        function print_slave(obj,~,evnt)
            postprocess = Postprocess_PhysicalProblem;
            res_file = evnt.AffectedObject.res_file;
            postprocess.print_slave(obj,res_file,obj.variables);
        end
    end
    
    methods (Abstract)
        preProcess(obj)
        computeVariables(obj)
        postProcess(obj)
    end
end
