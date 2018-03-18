classdef Physical_Problem < FEM
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        variables
        mesh
        dof
        problemID
        element
    end
    
    
    %% Restricted properties definition ===================================
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID,ptype)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID); 
            %% !! PREGUNTAR SI ELS HI SEMBLA BÉ A LA RESTA !!
            if nargin == 2
                obj.mesh.ptype = ptype;
            end
            obj.createGeometry(obj.mesh);
            obj.material = Material.create(obj.geometry,obj.mesh);
            obj.dof = DOF(problemID,obj.geometry,obj.mesh);
        end
        
        function preProcess(obj)
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.dof);
            obj.solver = Solver.create();
        end
        
        function computeVariables(obj)
            tol   = 1e-6;          
            for ifield = 1:obj.geometry(1).nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end

            transient = false;
              
             if transient
                  dt=0.01;
                  final_time = 1;
                  x = obj.solve_transient_problem(free_dof,tol,dt,final_time);  
              else
                  x = obj.solve_steady_problem(free_dof,tol);
              end
                       
             obj.variables = obj.element.computeVars(x);
        end
        
        function print(obj)
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,results);
        end
        
        function sol = solve_steady_problem(obj,free_dof,tol)

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
                try
            sol = x;
                catch
                    disp('eis')
                end
        end

        function sol = solve_transient_problem(obj,free_dof,tol,dt,final_time)
            total_free_dof= sum(free_dof);      
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
             x_n(:,istep)=x;
            end
            sol = x_n;
        end
        function postProcess(obj)
            %    ToDo
            % Inspire in TopOpt
            
        end
        
        function setDof(obj,dof)
            obj.dof = dof;
        end
        
        function setMatProps(obj,props)
            obj.element.material = obj.material.setProps(props);
        end
        
        function createGeometry(obj,mesh)
            
            if strcmp(mesh.ptype,'Stokes')
                obj.geometry=Geometry(mesh,'QUADRATIC');
                obj.geometry(2)=Geometry(mesh,'LINEAR');
                obj.geometry(1).nfields = 2;
            else
                obj.geometry=Geometry(mesh,'LINEAR');
            end
        end
    end   
end

