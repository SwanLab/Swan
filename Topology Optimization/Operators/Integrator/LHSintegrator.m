classdef LHSintegrator < handle
    
    properties (Access = private)
        quadrature
        interpolation        
        LHScells        
    end
    
    properties (Access = private)
        mesh
        globalConnec
        npnod
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
        end
        
        function LHS = compute(obj)
            lhs = obj.computeElementalLHS();
            LHS = obj.assembleMatrix(lhs);    
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.globalConnec  = cParams.globalConnec;
            obj.npnod         = cParams.npnod;
        end
        
       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature('LINEAR');            
           obj.quadrature = quad;
       end

        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;                        
        end        
        
        function lhs = computeElementalLHS(obj)
            shapes = obj.interpolation.shape;
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            lhs = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                dv(1,1,:) = dvolu(igaus,:);
                Ni = shapes(:,igaus);
                Nj = shapes(:,igaus);
                NiNj = Ni*Nj';
                Aij = bsxfun(@times,NiNj,dv);
                lhs = lhs + Aij;
            end
        end        
        
        function A = assembleMatrix(obj,aElem)
            connec = obj.globalConnec;
            ndofs  = obj.npnod;
            Ae     = aElem;
            nunkn1 = 1;
            nunkn2 = 1;
            nnode1 = size(connec,2);
            nnode2 = size(connec,2);
            A = sparse(ndofs,ndofs);
            for i = 1:nnode1*nunkn1
                nodeI = connec(:,i);                
                for j = 1:nnode2*nunkn2
                    nodeJ = connec(:,j);
                    a = squeeze(Ae(i,j,:));
                    A = A + sparse(nodeI,nodeJ,a,ndofs,ndofs);
                end
            end
        end        
        
    end
    
end