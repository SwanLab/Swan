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
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = public)
        material
        
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID);
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry = Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem,obj.mesh.connec,obj.geometry.cartd,obj.geometry.nnode,obj.mesh.coord);
            obj.dof = DOF(problemID,obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.mesh.scale);
        end
        
        function preProcess(obj)
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.bc,obj.dof,obj.dim);
            obj.solver  = Solver.create();
        end
        
        function computeVariables(obj)
            tol   = 1e-12;
            x     = zeros(length(obj.dof.free),1);
            inc_x = zeros(length(obj.dof.free),1);
            miter = 1e5;
            
            for incrm = 1:obj.element.nincr
                niter = 1;
                obj.element.cload = obj.element.cload + obj.element.fincr;
                [r,dr] = obj.element.computeResidual(inc_x);
                error = 1;
                while error > tol && niter <= miter

                    % Solve
                    inc_x = obj.solver.solve(dr,-r);
                    
                    % Updates
                    x = x + inc_x;
                    
                    % Compute new r & dr
                    [r,dr] = obj.element.computeResidual(inc_x);
                    
                    error = norm(r)/norm(obj.element.cload);
                    errcont(incrm,niter) = error;
                    niter = niter + 1;
                end
                nn(incrm) = niter-1;
                u = reshape(inc_x,2,[])';
                un(incrm) = u(1,1);
                
                xn(incrm) = obj.element.coord(4,1);
                fn(incrm) = obj.element.cload(obj.dof.free(1));

            end
            
%             % Convergence
%             Xcoord = log(errcont(end,1:niter-2));
%             Ycoord = log(errcont(end,2:niter-1));
%             [a,~] = polyfit(Xcoord,Ycoord,1);
%             fprintf('Convergence order, p: %d\nRatio, mu: %d\n\n',a(1),a(2));
%             
%             % Figures
%             figure;
%             [hAx,hLine1,hLine2] = plotyy(xn,fn,xn,nn,'plot','stairs');
%             hLine1.LineStyle = '-';
%             hLine1.Marker = '+';
%             xlabel('X displacement [m]')
%             ylabel(hAx(1),'Force [N]')
%             ylabel(hAx(2),'Iterations')
%             set(gcf,'Color','w')
%             set(hAx(1),'FontSize',15)
%             set(hAx(2),'FontSize',15)
            
            obj.variables = obj.element.computeVars(x);
        end
        
        function print(obj)
            iter = 1; % static
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,iter,results);
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
    %     methods (Access = protected, Static)s
    %     end
end

