classdef Physical_Problem < FEM
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        variables
        mesh
        dim
        dof
        bc
        problemID
        element
    end
    
    
    %% Restricted properties definition ===================================
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material

    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID); 
            obj.createGeometry(obj.mesh);
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.material = Material.create(obj.geometry,obj.mesh);
            obj.dof = DOF(problemID,obj.geometry,obj.dim,obj.mesh);%.scale,obj.nfields,obj.mesh.ptype,obj.interpolation_geometry,...
               % obj.mesh.nelem);
        end
        
        function preProcess(obj)
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.bc,obj.dof,obj.dim);
            obj.solver = Solver.create();
        end
        
        function computeVariables(obj)
            tol   = 1e-6;          
            for ifield = 1:obj.geometry.nfields
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
                while dot(r,r) > tol
                    inc_x = obj.solver.solve(dr,-r);
                    x = x0 + inc_x;
                    % Compute r 
                    r = obj.element.computeResidual(x,dr);
                    x0 = x;
                end
            sol=x;
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
        
        function [K, M] = computeKM(obj,job)
                        
            % !! Hyper-mega-ultra provisional !!
            
            dim_smooth.nnode = obj.geometry.nnode;
            dim_smooth.nunkn = 1;
            dim_smooth.nstre = 2;
            
            mesh_smooth = obj.mesh;
            mesh_smooth.ptype = 'DIFF-REACT';
            mesh_smooth.scale = 'MACRO';
            
            bc_smooth = obj.bc;
            bc_smooth.fixnodes = [];
            
            dof_smooth = DOF(obj.problemID,obj.geometry,obj.interpolation_variable,dim_smooth,mesh_smooth.scale,1,mesh_smooth.ptype,...
                obj.interpolation_geometry,obj.mesh.nelem);   
            dof_smooth.neumann = [];
            dof_smooth.dirichlet{1} = [];
            dof_smooth.neumann_values = [];
            dof_smooth.dirichlet_values{1} = [];
            dof_smooth.periodic_free = [];
            dof_smooth.periodic_constrained = [];
            dof_smooth.constrained = [];
            dof_smooth.free = setdiff(1:dof_smooth.ndof,dof_smooth.constrained);
            
            element_smooth = Element.create(mesh_smooth,obj.geometry,obj.material,obj.bc,dof_smooth,dim_smooth,1);
            
            [K] = element_smooth.computeStiffnessMatrix;
            [M] = element_smooth.computeMassMatrix(job);
        end
        function createGeometry(obj,mesh)
            obj.geometry=Geometry(mesh);
            if strcmp(mesh.ptype,'Stokes')
                obj.geometry(2)=Geometry(mesh);
                obj.geometry(1).interpolation.quadrature.computeQuadrature('QUADRATIC');
                obj.geometry(1).interpolation.isoparametric=obj.geometry(1).interpolation.set_isoparametric(mesh.geometryType,'QUADRATIC');
                obj.geometry(1).interpolation.compute_xpoints_T;
                obj.geometry(:).nfields = 2;
            end
        end
    end   
end

