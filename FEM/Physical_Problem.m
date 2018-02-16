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
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry = Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem,obj.mesh.connec,obj.geometry.cartd,obj.geometry.nnode,obj.mesh.coord);
            obj.bc = BC(obj.dim.nunkn,obj.problemID);
        end
        
        function preProcess(obj)
            % Create Objects
            obj.dof = DOF.create(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc,obj.mesh.scale);
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.bc,obj.dof,obj.dim);
%             obj.physicalVars = PhysicalVariables.create(obj.mesh.ptype,obj.mesh.pdim);
            obj.solver = Solver.create(obj.mesh.ptype);
        end
        
        function computeVariables(obj)
            tol   = 1e-6;
            x0 = zeros(length(obj.dof.vF),1);
            % Compute r & dr
            [r,dr] = obj.element.computeResidual(x0);
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r & dr
                [r,dr] = obj.element.computeResidual(x);
                x0 = x;
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
        
        function K = computeKsmooth(obj)
            
            dim_smooth.nnode=obj.geometry.nnode;
            dim_smooth.nunkn=1;
            dim_smooth.nstre=2;
            
            mesh_smooth = obj.mesh;
            mesh_smooth.ptype = 'THERMAL';
            mesh_smooth.scale = 'MACRO';
            
            bc_smooth = obj.bc;
            bc_smooth.fixnodes = [];
            
            dof_smooth = DOF.create(obj.geometry.nnode,obj.mesh.connec,dim_smooth.nunkn,obj.mesh.npnod,bc_smooth,mesh_smooth.scale);
            
            element_smooth = Element.create(mesh_smooth,obj.geometry,obj.material,obj.bc,dof_smooth,dim_smooth);

            [~,K] = element_smooth.computeResidual(zeros(dof_smooth.ndof,1));

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

