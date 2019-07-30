classdef Integrator_Simple < Integrator  
    
    properties (Access = private)
        globalConnec
        npnod
        Aelem
        Aglobal
        
        quadrature
        interpolation
        geometry
        
        innerShapes
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            %obj.globalConnec = obj.mesh.connec;
            obj.globalConnec = cParams.globalConnec;
            obj.npnod = cParams.npnod;
        end
        
        function A = computeLHS(obj)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.computeElementalMatrix();            
            obj.assembleMatrix();
            A = obj.Aglobal;
        end   
        
        function integrate()
            obj.initShapes();
            obj.innerShapes = obj.evaluateInnerShapes(F);
            obj.assembleLocally();
            A = obj.assembleIntegrand(obj.shapes);       
        end
        
    end
    
    methods (Access = private)
        
                function initShapes(obj)
            nelem = obj.backgroundMesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            obj.shapes = zeros(nelem,nnode);
                end
        
            function shapeValues = evaluateInnerShapes(obj,F1)
                interpolation = Interpolation.create(obj.backgroundMesh,'LINEAR');
                quadrature = obj.computeQuadrature(obj.backgroundMesh.geometryType);
                interpolation.computeShapeDeriv(quadrature.posgp);
                geometry = Geometry(obj.backgroundMesh,'LINEAR');
                geometry.computeGeometry(quadrature,interpolation);

                obj.innerShapes = zeros(size(obj.backgroundMesh.connec));
                for igauss = 1:quadrature.ngaus
                    obj.innerShapes = obj.innerShapes + interpolation.shape(:,igauss)'.*geometry.dvolu(:,igauss);
                end
            end
            
        function assembleLocally(obj)
            innerCells = obj.meshUnfitted.backgroundFullCells;
            obj.shapes(innerCells,:) = obj.innerShapes;
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
        
        function computeElementalMatrix(obj)            
            quad   = obj.quadrature;
            jacob  = obj.geometry.djacob;
            shapes = obj.interpolation.shape;
            ngaus  = quad.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            Ae = zeros(nnode,nnode,nelem);            
            for igaus = 1:ngaus
                for inode = 1:nnode
                    for jnode = 1:nnode
                        dJ = jacob(:,igaus);
                        w  = quad.weigp(igaus);
                        Ni = shapes(inode,igaus);
                        Nj = shapes(jnode,igaus);
                        Aij = squeeze(Ae(inode,jnode,:));
                        Ae(inode,jnode,:) = Aij + w*Ni*Nj*dJ;
                    end
                end
            end
            obj.Aelem = Ae;            
        end
    
        function assembleMatrix(obj)
            connec = obj.globalConnec;
            ndofs  = obj.npnod;
            Ae     = obj.Aelem;
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
            obj.Aglobal = A;
        end
        
    end
    
end