classdef Geometry_Surface < Geometry
    
    properties (Access = public)
        normalVector        
    end
    
    properties (Access = private)
       drDtxi 
       jacobian
    end
    
    methods (Access = public)
        
        function obj = Geometry_Surface(cParams)
            obj.permutation = [3 2 1];                        
            obj.init(cParams);
        end
        
        function computeGeometry(obj,quad,interpV)
            obj.initGeometry(interpV,quad);
            obj.computeDrDtxi();
            obj.computeNormals();
            obj.computeJacobian();
            obj.computeDvolu();
        end
        
    end
    
    methods (Access = private)

        function computeDrDtxi(obj)
            nEmb    = 2;
            nDime   = obj.mesh.ndim;
            nGaus   = obj.quadrature.ngaus;
            nElem   = obj.mesh.nelem;
            nNode   = obj.mesh.nnode;            
            xp      = obj.coordElem;
            deriv   = obj.mesh.interpolation.deriv(:,:,:);
            dShapes = permute(deriv,[1 3 2]);
            obj.drDtxi = zeros(nGaus,nElem,nEmb,nDime);            
            for idime = 1:nDime
                dxDtxi = zeros(nElem,nEmb,nGaus);
                for inode = 1:nNode
                    dShape(1,:,:) = dShapes(:,:,inode);
                    xV(:,1)= xp(:,inode,idime);
                    dxDtxi = dxDtxi + bsxfun(@times,xV,dShape);
                end
                obj.drDtxi(:,:,:,idime) = permute(dxDtxi,[3 1 2]);
            end            
        end
        
        function computeNormals(obj)
            nDime   = obj.mesh.ndim;
            nGaus   = obj.quadrature.ngaus;
            nElem   = obj.mesh.nelem;            
            obj.normalVector = zeros(nGaus,nElem,nDime);            
            DxDtxi = obj.drDtxi(:,:,1,1);
            DxDeta = obj.drDtxi(:,:,2,1);
            DyDtxi = obj.drDtxi(:,:,1,2);
            DyDeta = obj.drDtxi(:,:,2,2);
            DzDtxi = obj.drDtxi(:,:,1,3);
            DzDeta = obj.drDtxi(:,:,2,3);            
            obj.normalVector(:,:,1) = DyDtxi.*DzDeta - DzDtxi.*DyDeta;
            obj.normalVector(:,:,2) = DzDtxi.*DxDeta - DxDtxi.*DzDeta;
            obj.normalVector(:,:,3) = DxDtxi.*DyDeta - DyDtxi.*DxDeta;            
        end
        
        function computeJacobian(obj)
            nx = obj.normalVector(:,:,1);
            ny = obj.normalVector(:,:,2);
            nz = obj.normalVector(:,:,3);
            jac = sqrt(nx.^2 + ny.^2 + nz.^2);
            obj.jacobian = jac;
        end
        
        function computeDvolu(obj)
            jac = obj.jacobian;
            w(:,1) = obj.quadrature.weigp;
            dv =  bsxfun(@times,w,jac);
            obj.dvolu = dv';
        end
        
    end
    
end