classdef Element_DiffReact < Element
    %Element_DiffReact Summary of this class goes here
    %   Detailed explanation goes here
    
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
        mesh
        K
        M
        Mr
        epsilon
        interpolation_u
    end
    
    properties (Access = private)
        nstre
        addRobinTerm
        boundaryMesh
    end
    
    methods %(Access = ?Physical_Problem)
        function obj = Element_DiffReact(mesh,geometry,material,dof,scale,addRobinTerm,bcType,interp,boundaryMesh)
            obj.mesh = mesh;
            obj.addRobinTerm = addRobinTerm;
            obj.bcType = bcType;
            obj.initElement(geometry,mesh,material,dof,scale,interp);
            obj.nstre = 2;
            obj.nfields = 1;
            obj.interpolation_u = interp{1};
            obj.boundaryMesh = boundaryMesh;
            obj.computeStiffnessMatrix();
            obj.computeMassMatrix();
            obj.computeBoundaryMassMatrix();
        end
        
        function obj = setEpsilon(obj,epsilon)
            obj.epsilon = epsilon;
        end
        
        function LHS = computeLHS(obj)
            if obj.addRobinTerm
                LHS = obj.epsilon^2*obj.K + obj.M + (obj.epsilon)*obj.Mr;
            else
                LHS = obj.epsilon^2*obj.K + obj.M;
                LHS = obj.bcApplier.fullToReducedMatrix(LHS);
            end
        end
        

    end
    
    methods (Access = private)
        
        function computeStiffnessMatrix(obj)
            Ke = obj.computeElementalStiffnessMatrix();
            Kg = obj.AssembleMatrix(Ke,1,1); % !!
            obj.K = Kg;
        end
        
        function computeMassMatrix(obj)
            Me = obj.computeElementalMassMatrix();
            Mg = obj.AssembleMatrix(Me,1,1); % !!
            obj.M = Mg;
        end
        
        function computeBoundaryMassMatrix(obj)
            if obj.addRobinTerm
                cParams = obj.createIntegratorParams();
                integrator = Integrator.create(cParams);
                obj.Mr = integrator.computeLHS();
            end
        end
        
        function Ke = computeElementalStiffnessMatrix(obj)
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
            ndof  = obj.dof.nunkn*obj.nnode;
            ngaus = obj.quadrature.ngaus;
            Ke = zeros(ndof,ndof,obj.nelem);
            for igaus = 1:ngaus
                dShapeDx = obj.geometry.cartd(:,:,:,igaus);
                Bmat     = obj.computeB(obj.dof.nunkn,obj.nelem,obj.nnode,dShapeDx);
                for istre = 1:obj.nstre
                    BmatI = Bmat(istre,:,:);
                    BmatJ = permute(Bmat(istre,:,:),[2 1 3]);
                    dNdN = bsxfun(@times,BmatJ,BmatI);
                    dv(1,1,:) = obj.geometry.dvolu(:,igaus);
                    inc = bsxfun(@times,dv,dNdN);
                    Ke = Ke + inc;
                end
            end            
        end
        
        function Me = computeElementalMassMatrix(obj)
            obj.quadrature.computeQuadrature('QUADRATICMASS');
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);            
            shapes = obj.interpolation_u.shape;
            dvolu  = obj.geometry.dvolu;
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            Me = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                dv(1,1,:) = dvolu(:,igaus);
                Ni = shapes(:,igaus);
                Nj = shapes(:,igaus);
                NiNj = Ni*Nj';
                Mij = bsxfun(@times,NiNj,dv);
                Me = Me + Mij;
            end            
            obj.quadrature.computeQuadrature('LINEAR');
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation_u);
        end        
        
        function params = createIntegratorParams(obj)
            params.type = 'COMPOSITE';
            params.npnod = obj.mesh.npnod;
            
            s.backgroundMesh = obj.mesh;
            s.dimension = 1:s.backgroundMesh.ndim;
            s.type = 'FromReactangularBox';
            bC = BoundaryMeshCreator.create(s);
            bMeshes = bC.create();
            
            bMeshes = obj.boundaryMesh;
            
            nBoxFaces = numel(bMeshes);
            
            
            for iMesh = 1:nBoxFaces
                boxFaceMesh = bMeshes{iMesh};
                cParams.mesh = boxFaceMesh.mesh;
                cParams.type = 'SIMPLE';
                cParams.globalConnec = boxFaceMesh.globalConnec;
                cParams.npnod        = obj.mesh.npnod;
                cParams.geometryType = obj.mesh.type;
                params.compositeParams{iMesh} = cParams;
            end
        end
        
    end
    
    methods(Access = protected) % Only the child sees the function
        function FextSuperficial = computeSuperficialFext(obj)
            FextSuperficial = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
        
        function FextVolumetric = computeVolumetricFext(obj)
            FextVolumetric = zeros(obj.nnode*obj.dof.nunkn,1,obj.nelem);
        end
    end
    
    methods (Static)
        function [B] = computeB(nunkn,nelem,nnode,cartd)
            B = zeros(2,nnode*nunkn,nelem);
            for inode=1:nnode
                j = nunkn*(inode-1)+1;
                B(1,j,:)=cartd(1,inode,:);
                B(2,j,:)=cartd(2,inode,:);
            end
        end
    end
end


