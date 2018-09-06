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
            mesh = Mesh_GiD(problemID); % Mesh defined twice, but almost free
            switch mesh.ptype
                case 'ELASTIC'
                    switch mesh.scale                        
                        case 'MACRO'
                            obj = Elastic_Problem(problemID);
                        case 'MICRO'
                            obj = Elastic_Problem_Micro(problemID);
                    end                    
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
        function sol = solve_steady_nonlinear_problem(obj,free_dof,tol)
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
        
        function sol = solve_transient_nonlinear_problem(obj,free_dof,tol,dt,final_time)
            total_free_dof = sum(free_dof);
            x_n(:,1) = zeros(total_free_dof,1);
            x0 = zeros(total_free_dof,1);
            
            dr = obj.element.computedr(dt);
            
            for istep = 2: final_time/dt
                u_previous_step = x_n(1:free_dof(1),istep-1);
                
                r = obj.element.computeResidual(x0,dr,u_previous_step);
                while dot(r,r) > tol
                    inc_x = obj.solver.solve(dr,-r);
                    x = x0 + inc_x;
                    % Compute r
                    r = obj.element.computeResidual(x,dr,u_previous_step);
                    x0 = x;
                end
                x_n(:,istep) = x;
            end
            sol = x_n;
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
