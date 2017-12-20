classdef Physical_Problem < FEM
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)
        variables
        mesh
        dim
        problemID
    end
    
    %% Private properties definition ======================================
    properties (Access = protected)
        bc
        material
        element
        dof
        physicalVars
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID);
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry=Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem); 
        end
        
        function preProcess(obj)
            % Create Objects
            obj.bc = BC(obj.dim.nunkn,obj.problemID);
            obj.dof = DOF(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.fixnodes);
            obj.element = Element.create(obj.mesh.ptype,obj.mesh.pdim);
            obj.physicalVars = PhysicalVariables.create(obj.mesh.ptype,obj.mesh.pdim);
            obj.solver = Solver_Dirichlet_Conditions;
        end
        
        function computeVariables(obj)
            obj.element.computeLHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry,obj.material);
            obj.element.computeRHS(obj.dim.nunkn,obj.mesh.nelem,obj.geometry.nnode,obj.bc,obj.dof.idx);
            
            % Assembly
            [obj.LHS,obj.RHS] = obj.Assemble(obj.element,obj.geometry.nnode,obj.dim.nunkn,obj.dof);
            
            % Solver
            sol = obj.solver.solve(obj.LHS,obj.RHS,obj.dof,obj.bc.fixnodes);
            obj.variables = obj.physicalVars.computeVars(sol,obj.dim,obj.geometry,obj.mesh.nelem,obj.dof.idx,obj.element,obj.material);
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
        function Msmooth=computeMass(obj,job)
            meshMass=obj.mesh;
            meshMass.geometryType='Triangle_Linear_Mass';
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
        function StifMat=computeKsmooth(obj)
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
    methods (Access = protected, Static)
        function [LHS,RHS] = Assemble(element,nnode,nunkn,dof)            
            % Compute LHS
            LHS = sparse(dof.ndof,dof.ndof);
            for i = 1:nnode*nunkn
                for j = 1:nnode*nunkn
                    a = squeeze(element.LHS(i,j,:));
                    LHS = LHS + sparse(dof.idx(i,:),dof.idx(j,:),a,dof.ndof,dof.ndof);
                end
            end            
            % Compute RHS
            RHS = zeros(dof.ndof,1);
            for i = 1:length(dof.idx(:,1)) % nnode*nunkn
                b = squeeze(element.RHS(i,1,:));
                ind = dof.idx(i,:);
                RHS(ind) = b;
            end
        end
    end
end

