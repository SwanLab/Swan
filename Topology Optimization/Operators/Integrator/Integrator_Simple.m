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
        innerToBackground
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.globalConnec = cParams.globalConnec;
            obj.innerToBackground = cParams.innerToBackground;
            obj.npnod = cParams.npnod;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
        end
        
        function LHS = computeLHS(obj)
            obj.computeElementalLHS();
            obj.assembleMatrix();
            LHS = obj.LHS;
        end
        
        function RHS = integrate(obj,F)
            obj.initShapes();
            obj.computeElementalRHS(F);
            obj.assembleSubcellsInCells();
            RHS = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        function initShapes(obj)
            nelem = obj.backgroundMesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            obj.RHScells = zeros(nelem,nnode);
        end
        
        function computeElementalRHS(obj,F1)
            jacob  = obj.geometry.djacob;
            shapes = obj.interpolation.shape;
            weight = obj.quadrature.weigp;
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            f      = zeros(nelem,nnode);
            for igaus = 1:ngaus
                for iNode = 1:nnode
                    dJ = jacob(:,igaus);
                    w = weight(igaus);
                    Ni = shapes(iNode,igaus);
                    f(:,iNode) = f(:,iNode) + w*Ni*dJ;
                end
            end
            obj.RHSsubcells = f;
        end
        
        function assembleSubcellsInCells(obj)
            innerCells = obj.innerToBackground;
            obj.RHScells(innerCells,:) = obj.RHSsubcells;
        end
        
        function f = assembleIntegrand(obj)
            integrand = obj.RHScells;
            npnod  = obj.backgroundMesh.npnod;
            nnode  = obj.backgroundMesh.nnode;
            connec = obj.backgroundMesh.connec;
            f = zeros(npnod,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[npnod,1],@sum,0);
            end
        end
        
        function createQuadrature(obj)
            quad = obj.computeQuadrature(obj.mesh.geometryType);
            obj.quadrature = quad;
        end
        
        function createInterpolation(obj)
            quad = obj.quadrature;
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(quad.posgp);
            obj.interpolation = int;
        end
        
        function createGeometry(obj)
            quad = obj.quadrature;
            int  = obj.interpolation;
            geom = Geometry(obj.mesh,'LINEAR');
            geom.computeGeometry(quad,int);
            obj.geometry = geom;
        end
        
        function computeElementalLHS(obj)
            jacob  = obj.geometry.djacob;
            shapes = obj.interpolation.shape;
            weight  = obj.quadrature.weigp;
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            Ae = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                for inode = 1:nnode
                    for jnode = 1:nnode
                        dJ = jacob(:,igaus);
                        w  = weight(igaus);
                        Ni = shapes(inode,igaus);
                        Nj = shapes(jnode,igaus);
                        Aij = squeeze(Ae(inode,jnode,:));
                        Ae(inode,jnode,:) = Aij + w*Ni*Nj*dJ;
                    end
                end
            end
            obj.LHScells = Ae;
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