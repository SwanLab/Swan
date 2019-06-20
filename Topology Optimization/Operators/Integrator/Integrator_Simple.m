classdef Integrator_Simple < Integrator  
    
    properties (Access = private)
        globalConnec
        npnod
        Aelem
        Aglobal
    end
    
    methods (Access = public)
        
        function obj = Integrator_Simple(cParams)
            obj.init(cParams)
            obj.globalConnec = cParams.globalConnec;
            obj.npnod = cParams.npnod;
        end
        
        function A = computeLHS(obj)
            obj.computeElementalMatrix();            
            obj.assambleMatrix();
            A = obj.Aglobal;
        end                
        
    end
    
    methods (Access = private)
        
        function computeElementalMatrix(obj)
            interpolation = Interpolation.create(obj.mesh,'LINEAR');
            quadrature = obj.computeQuadrature(obj.mesh.geometryType);
            interpolation.computeShapeDeriv(quadrature.posgp);
            geometry = Geometry(obj.mesh,'LINEAR');
            geometry.computeGeometry(quadrature,interpolation);
            nelem = obj.mesh.nelem;
            Ae = zeros(interpolation.nnode,interpolation.nnode,nelem);            
            for igaus = 1:quadrature.ngaus
                for inode = 1:interpolation.nnode
                    for jnode = 1:interpolation.nnode
                        Ae(inode,jnode,:) = squeeze(Ae(inode,jnode,:)) + quadrature.weigp(igaus)*interpolation.shape(inode,igaus)...
                            *interpolation.shape(jnode,igaus)*geometry.djacob(:,igaus);
                    end
                end
            end
            obj.Aelem = Ae;            
        end
        
        function assambleMatrix(obj)
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