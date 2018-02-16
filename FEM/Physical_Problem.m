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
        
    end
    
    
    %% Restricted properties definition ===================================
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = public)
        material
        element
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID);
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry = Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem,obj.mesh.connec,obj.geometry.cartd,obj.geometry.nnode,obj.mesh.coord);
            obj.bc = BC(obj.dim.nunkn,obj.problemID);
        end
        
        function preProcess(obj)
            % Create Objects
            obj.dof = DOF(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.fixnodes);
            obj.element = Element.create(obj.mesh,obj.dim,obj.geometry,obj.material,obj.bc,obj.dof);
            obj.solver = Solver.create(obj.mesh.ptype);
        end
        
        function computeVariables(obj)
            tol   = 1e-12;
            x     = zeros(length(obj.dof.vL),1);
            r     = 0;
            resid = 0;
            cload = 0;

            miter = 1e4;
            nincr = 100;
            
            % Compute r & dr
            fincr = 0.5*obj.element.Fext/nincr;
            for incrm = 1:nincr
                niter = 1;
                
                % Set current load & residual
                cload = cload + fincr;
                resid = resid - fincr;
                r = r - fincr(obj.dof.vL);
                
                while dot(r,r) > tol && niter <= miter
                    [~,dr,K] = obj.element.computeResidual(x);
                    inc_x = obj.solver.solve(dr,-r);
                    
                    x = x + inc_x;
                    
                    % Updates
                    obj.element.updateCoord(inc_x);
                    obj.element.updateCartd(obj.mesh.pdim);
                    
                    % Compute fint
                    fint = obj.element.computeInternal();
                    
                    % Compute residual
                    r = fint - cload;
                    r = r(obj.dof.vL);
                    niter = niter + 1;
                end
            end
            obj.variables = obj.element.computeVars(x);
        end
        
        function postProcess(obj)
            iter = 1; % static
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,iter,results);
        end
        
        function setMatProps(obj,props)
            obj.material = obj.material.setProps(props);
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
            lnods=obj.mesh.connec';
            emat = zeros(geom.nnode,geom.nnode,obj.mesh.nelem);
            for igaus=1:geom.ngaus
                for inode=1:geom.nnode
                    for jnode=1:geom.nnode
                        emat(inode,jnode,:)=squeeze(emat(inode,jnode,:)) + geom.weigp(igaus)*geom.shape(inode,igaus)*geom.shape(jnode,igaus)*geom.djacob(:,igaus);
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
                    Msmooth = Msmooth + sparse(lnods(inode,:),1,elumped(inode,:),npnod,1);
                end
            elseif (job==2)
                
                Msmooth = sparse(obj.mesh.npnod,obj.mesh.npnod);
                for k=1:geom.nnode
                    for l=1:geom.nnode
                        vmass = squeeze(emat(k,l,:));
                        Msmooth = Msmooth + sparse(lnods(k,:),lnods(l,:),vmass,obj.mesh.npnod,obj.mesh.npnod);
                    end
                end
                
            end
        end
        
        function StifMat = computeKsmooth(obj)
            StifMat=sparse(obj.mesh.npnod,obj.mesh.npnod);
            nnode=obj.geometry.nnode;
            nunkn=1;
            nstre=2;
            element_smooth=Element.create('THERMAL',obj.mesh.pdim);
            element_smooth.computeLHS(nunkn,nstre,obj.mesh.nelem,obj.geometry,obj.material);
            estiff=element_smooth.LHS;
            lnods=obj.mesh.connec';
            for a=1:nnode
                for i=1:nunkn
                    idx(nunkn*a-nunkn+i,:) = nunkn.*lnods(a,:)-nunkn+i;
                end
            end
            for k=1:nnode*nunkn
                for l=1:nnode*nunkn
                    vestiff = squeeze(estiff(k,l,:));
                    StifMat = StifMat + sparse(idx(k,:),idx(l,:),vestiff,obj.mesh.npnod,obj.mesh.npnod);
                end
            end
        end
    end
    
    %% Private methods definition =========================================
    %     methods (Access = protected, Static)s
    %     end
end

