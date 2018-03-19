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
            x     = zeros(length(obj.dof.free),1); % acumulada
            x_i   = zeros(length(obj.dof.free),1);
            miter = 1e5;
            lambda_i = 0;
            
            
            for incrm = 1:obj.element.nincr     % i
                niter = 1;                      % k

                lambda_i = lambda_i + 1/obj.element.nincr;
                
                [r,dr,K,fint] = obj.element.computeResidual(x_i,lambda_i);
                error = 1;
                lambda_k = lambda_i;
                x_k = x_i;
                while error > tol && niter <= miter
                    % *****************************************************
                    
                    R = fint(obj.dof.constrained) - lambda_k*obj.element.fext(obj.dof.constrained);
                    uR = K(obj.dof.constrained,obj.dof.constrained)\R;
                    uF = -K(obj.dof.free,obj.dof.free)\obj.element.fext(obj.dof.free);
                    
                    inc_lambda_k = lambda_k - lambda_i; % inc_lambda_k=0 en primera iter de k
                    inc_xk = x_k - x_i;
                    
                    Psi = 1;
                    
                    s_k = [inc_xk;inc_lambda_k*Psi*obj.element.fext];
                    s = norm(s_k);

                    a1 = uF'*uF + Psi^2*obj.element.fext'*obj.element.fext;
                    a2 = 2*uF'*(inc_xk+uR) + 2*inc_lambda_k*Psi^2*obj.element.fext'*obj.element.fext;
                    a3 = uR'*(2*inc_xk+uR) + inc_xk'*inc_xk - s^2;

                    poly2 = [a1 a2 a3];
                    gamma = roots(poly2);
                    
                    u_1 = uR + gamma(1)*uF;
                    u_2 = uR + gamma(2)*uF;
                    
                    s_k1_1 = [inc_xk+u_1; (inc_lambda_k+gamma(1))*Psi*obj.element.fext];
                    s_k1_2 = [inc_xk+u_2; (inc_lambda_k+gamma(2))*Psi*obj.element.fext];
                    
                    theta_1= acosd(dot(s_k,s_k1_1)/s^2);
                    theta_2= acosd(dot(s_k,s_k1_2)/s^2);
                    theta = [theta_1 theta_2];
                    
                    [~,j] = min(theta);
                    gamma = gamma(j);
                    
                    u_k = uR + gamma*uF;
                    lambda_k = lambda_k + gamma;
                    inc_lambda_k = inc_lambda_k + gamma;
                    inc_xk = inc_xk + u_k;
                    x_k = x_k + u_k;
                    
                    % *****************************************************
                    % Solve
%                     x_k = obj.solver.solve(dr,-r);
                    
                    % Updates
                    x = x + x_k;
                    
                    % Compute new r & dr
                    [r,dr,K,fint] = obj.element.computeResidual(x_k,lambda_k);
                    
                    error = norm(r)/norm(obj.element.fext);
                    errcont(incrm,niter) = error;
                    niter = niter + 1;
                end
%                 x_i = x_k;
%                 lambda_i = lambda_k;
                
                
                global test
                if exist('test','var') == 1
                    nn(incrm) = niter-1;

                    gid_elem_pos = 4; %6
                    gid_dime_pos = 1; % X Y Z 2

                    xn(incrm) = obj.element.coord(gid_elem_pos,gid_dime_pos);
                    gid_id = gid_elem_pos*obj.geometry.ndime-obj.geometry.ndime+gid_dime_pos;
                    fn(incrm) = obj.element.fext(gid_id);

                    %% Plot
                    xlabel('X displacement [m]')
                    set(gcf,'Color','w')

                    % Left
                    yyaxis left
                    plot(xn,fn,'*')
                    ylabel('Force [N]')
                    set(gca,'FontSize',15)

                    % Right
                    yyaxis right
                    plot(xn,nn,'o-')
                    ylabel('Iterations')
                    set(gca,'FontSize',15)

                    % Update
                    drawnow;
                end
            end
            
            % Convergence
            if exist('test','var') == 1
            Xcoord = log(errcont(end,1:niter-2));
            Ycoord = log(errcont(end,2:niter-1));
            [a,~] = polyfit(Xcoord,Ycoord,1);
            fprintf('Convergence order, p: %d\nRatio, mu: %d\n\n',a(1),a(2));
            end
                        
            obj.variables = obj.element.computeVars(x);
        end
        
        function print(obj)
            for iter = 1:obj.element.nincr
                postprocess = Postprocess_PhysicalProblem();
                results.physicalVars = obj.variables;
                postprocess.print(obj,obj.problemID,iter,results);
            end
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

