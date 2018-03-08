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
        interpolation_geometry
        interpolation_variable
        quadrature
        geometry_variable
        nfields
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
            obj.interpolation_geometry = Interpolation.create('mesh');
            obj.interpolation_geometry.compute(obj.mesh);
            
            
            if strcmp(obj.mesh.ptype,'Stokes') == 1 
                order= {'QUADRATIC','LINEAR'};
                obj.nfields = 2;
            else
                order= {obj.interpolation_geometry.order};
                obj.nfields = 1;
            end

            obj.quadrature = Quadrature (obj.mesh.geometryType,strjoin(order(1)));
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
                
                for i=1:obj.nfields
                interpolation_variable(i) = Interpolation.create ('variable');
                interpolation_variable(i).compute(obj.interpolation_geometry,strjoin(order(i)));
                geometry_variable(i) = Geometry(interpolation_variable(i),obj.quadrature,obj.mesh.nelem); 
%                 dof(i) = DOF(problemID,geometry_variable(i).nnode,interpolation_variable(i).T,obj.dim.nunkn(i),interpolation_variable(i).npnod,obj.mesh.scale);
                end
            obj.interpolation_variable = interpolation_variable;
            obj.geometry_variable = geometry_variable;
%             obj.dof=dof;

%             obj.quadrature = Quadrature (obj.mesh.geometryType,strjoin(order(1)));
%             obj.interpolation_variable=Interpolation.create ('variable');
%             obj.interpolation_variable.compute(obj.interpolation_geometry,order,obj.nfields);
%             obj.geometry = Geometry(obj.interpolation_variable,obj.quadrature,obj.mesh.nelem,obj.nfields); 

            
%             obj.interpolation_variable = interpolation_variable; 
            
%             obj.geometry = geometry;
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem,obj.mesh.connec,obj.geometry_variable(1).cartd,obj.geometry_variable(1).nnode,obj.mesh.coord);
            obj.dof = DOF(problemID,obj.geometry_variable,obj.interpolation_variable,obj.dim,obj.mesh.scale,obj.nfields,obj.mesh.ptype,obj.interpolation_geometry,obj.mesh.nelem);
            
        end
        
        function preProcess(obj)
            obj.element = Element.create(obj.mesh,obj.geometry_variable,obj.material,obj.bc,obj.dof,obj.dim,obj.nfields);
            obj.solver = Solver.create();
        end
        
        function computeVariables(obj)
            tol   = 1e-6;
          
            for ifield = 1:obj.nfields
                free_dof(ifield) = length(obj.dof.free{ifield});
            end

            transient = true;
              
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
        
        function Msmooth = computeMass(obj,job)
            meshMass=obj.mesh;
            switch obj.geometry.type
                case 'TRIANGLE'
                    meshMass.geometryType='Triangle_Linear_Mass';
                case 'QUADRILATERAL'
                    meshMass.geometryType='Quad_Mass';
            end
            geom=Geometry(meshMass);
            dirichlet_data=obj.mesh.connec';
            emat = zeros(geom.nnode,geom.nnode,obj.mesh.nelem);
            for igaus=1:geom.ngaus
                for inode=1:geom.nnode
                    for jnode=1:geom.nnode
                        emat(inode,jnode,:)=squeeze(emat(inode,jnode,:)) + geom.weigp(igaus)*geom.shape(inode,igaus)*geom.shape(jnode,igaus)*geom.djacb(:,igaus);
                    end
                end
            end
            
            if (job==1)
                % lumped mass matrix
                elumped = zeros(geom.nnode,obj.mesh.nelem);
                Msmooth = zeros(geom.nnode,1);
                [nproc,coeff] = nprocedure(etype,nnode);
                if (nproc==1)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
                        end
                    end
                elseif (nproc==2)
                    for inode=1:nnode
                        for jnode=1:nnode
                            elumped(inode,:)=elumped(inode,:)+squeeze(emat(inode,jnode,:))';
                        end
                        elumped(inode,:)=elumped(inode,:)*coeff(inode);
                    end
                end
                for inode=1:nnode
                    Msmooth = Msmooth + sparse(dirichlet_data(inode,:),1,elumped(inode,:),npnod,1);
                end
            elseif (job==2)
                
                Msmooth = sparse(obj.mesh.npnod,obj.mesh.npnod);
                for k=1:geom.nnode
                    for l=1:geom.nnode
                        vmass = squeeze(emat(k,l,:));
                        Msmooth = Msmooth + sparse(dirichlet_data(k,:),dirichlet_data(l,:),vmass,obj.mesh.npnod,obj.mesh.npnod);
                    end
                end
                
            end
        end
        
        function K = computeKsmooth(obj)
            
            % Hyper-mega-ultra provisional
            
            dim_smooth.nnode=obj.geometry.nnode;
            dim_smooth.nunkn=1;
            dim_smooth.nstre=2;
            
            mesh_smooth = obj.mesh;
            mesh_smooth.ptype = 'THERMAL';
            mesh_smooth.scale = 'MACRO';
            
            bc_smooth = obj.bc;
            bc_smooth.fixnodes = [];
            
            dof_smooth = DOF(obj.problemID,obj.geometry.nnode,obj.mesh.connec,...
                             dim_smooth.nunkn,obj.mesh.npnod,mesh_smooth.scale);
           
            dof_smooth.neumann = [];
            dof_smooth.dirichlet = [];
            dof_smooth.neumann_values = [];
            dof_smooth.dirichlet_values = [];
            dof_smooth.periodic_free = [];
            dof_smooth.periodic_constrained = [];
            dof_smooth.constrained = [];
            dof_smooth.free = setdiff(1:dof_smooth.ndof,dof_smooth.constrained);
            
            
            element_smooth = Element.create(mesh_smooth,obj.geometry,obj.material,obj.bc,dof_smooth,dim_smooth);

            [K] = element_smooth.computeStiffnessMatrix;

        end
    end
    
    %% Private methods definition =========================================
    methods (Access = protected, Static)
        %         function [LHS,RHS] = Assemble(element,nnode,nunkn,dof,bc)
        %             % Compute LHS
        %             LHS = sparse(dof.ndof,dof.ndof);
        %             for i = 1:nnode*nunkn
        %                 for j = 1:nnode*nunkn
        %                     a = squeeze(element.LHS(i,j,:));
        %                     LHS = LHS + sparse(dof.idx(i,:),dof.idx(j,:),a,dof.ndof,dof.ndof);
        %                 end
        %             end
        %             LHS = 1/2 * (LHS + LHS');
        %             % Compute RHS
        %             RHS = zeros(dof.ndof,1);
        %             for i = 1:nnode*nunkn
        %                 b = squeeze(element.Fext(i,1,:));
        %                 ind = dof.idx(i,:);
        %                 RHS = RHS + sparse(ind,1,b',dof.ndof,1);
        %             end
        %
        %             %Compute Global Puntual Forces
        %             if ~isempty(bc.iN)
        %                 FextPoint = zeros(dof.ndof,1);
        %                 FextPoint(bc.iN) = bc.neunodes(:,3);
        %                 RHS = RHS + FextPoint;
        %             end
        %         end
    end
end

