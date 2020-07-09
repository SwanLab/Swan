classdef Integrator_Simple < Integrator
    
    properties (Access = private)
        globalConnec
        backgroundMesh
        npnod
        LHScells
        LHS
        
        quadrature
        interpolation
        geometry
        
        RHScells
        RHSsubcells
        Fnodal
        xGauss
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.npnod          = cParams.npnod;            
            obj.globalConnec   = cParams.globalConnec;
            
            
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end
        
        function LHS = computeLHS(obj)
            obj.computeElementalLHS();
            obj.assembleMatrix();
            LHS = obj.LHS;
        end
        
        function norm = computeL2norm(obj,f)
            fv(:,1) = f(:);
            M = obj.LHS;
            norm = fv'*M*fv;
        end
        
        function RHS = integrate(obj,F)
            obj.Fnodal = F;
            obj.initShapes();
            obj.computeElementalRHS();
            RHS = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        function initShapes(obj)
            nelem = obj.backgroundMesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            obj.RHScells = zeros(nelem,nnode);
        end
        
        function createQuadrature(obj)
            quad = obj.computeQuadrature(obj.mesh.geometryType);
            obj.quadrature = quad;
        end
        
        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            obj.interpolation = int;
        end
        
        function createGeometry(obj)
            quad = obj.quadrature;
            int  = obj.interpolation;
            s.mesh = obj.mesh;
            geom = Geometry.create(s);
            geom.computeGeometry(quad,int);
            obj.geometry = geom;
        end
        
        function computeElementalRHS(obj,F1)
            shapes = obj.interpolation.shape;
            dvolu  = obj.geometry.dvolu;
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            f      = zeros(nelem,nnode);
            for igaus = 1:ngaus
                dv = dvolu(:,igaus);
                Ni(1,:) = shapes(:,igaus);
                inc = bsxfun(@times,dv,Ni);
                f = f + inc;
            end
            obj.RHScells = f;
        end
        
%         function computeElementalRHS(obj)
%             obj.computeUnfittedGaussPoints();            
%             %obj.computeShapeFunctions();
%             int = obj.integrateFwithShapeFunction();
%             obj.RHScells = int;
%         end
        
        function computeUnfittedGaussPoints(obj)
            q = obj.quadrature;
            m = obj.mesh;
            obj.xGauss = m.computeXgauss(q.posgp);
        end
        
        function shapes = computeShapeFunctions(obj)
            m = obj.backgroundMesh;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);            
        end          
        
        function int = integrateFwithShapeFunction(obj)
            Fgauss  = obj.computeFinGaussPoints();
            dV      = obj.computeDvolume;
            fdV     = (Fgauss.*dV);
            shapes  = obj.computeShapeFunctions();
            int = obj.initIntegrand();
            for igaus = 1:obj.quadrature.ngaus
                fdv   = fdV(igaus,:);
                shape = shapes(:,:,igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            int = transpose(int);
        end   
        
        function int = initIntegrand(obj)
            nelem = obj.mesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            int = zeros(nnode,nelem);            
        end        
        
        function dV = computeDvolume(obj)
            q  = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);            
        end        
        
        function fG = computeFinGaussPoints(obj)
            f = createFeFunction(obj);
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);            
        end
        
        function f = createFeFunction(obj)            
            m = obj.mesh;
            s.fNodes = obj.Fnodal;
            s.mesh = m;
            f = FeFunction(s);
        end           
        
        function f = assembleIntegrand(obj)
            integrand = obj.RHScells;
            ndofs  = obj.backgroundMesh.npnod;
            connec = obj.globalConnec;
            nnode  = size(connec,2);
            f = zeros(ndofs,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[ndofs,1],@sum,0);
            end
        end
        
        function computeElementalLHS(obj)
            shapes = obj.interpolation.shape;
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
                Aij = bsxfun(@times,NiNj,dv);
                Me = Me + Aij;
            end
            obj.LHScells = Me;
        end
        
        function assembleMatrix(obj)
            connec = obj.globalConnec;
            ndofs  = obj.npnod;
            Ae     = obj.LHScells;
            nunkn1 = 1;
            nunkn2 = 1;
            nnode1 = size(connec,2);
            nnode2 = size(connec,2);
            idx1 = connec';
            idx2 = connec';
            tic
            A = sparse(ndofs,ndofs);
            for i = 1:nnode1*nunkn1
                for j = 1:nnode2*nunkn2
                    a = squeeze(Ae(i,j,:));
                    A = A + sparse(idx1(i,:),idx2(j,:),a,ndofs,ndofs);
                end
            end
            obj.LHS = A;
        end
        
    end
    
end