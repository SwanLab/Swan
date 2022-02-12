classdef LHSintegrator < handle

    properties (Access = protected)
        mesh   
        quadrature
        interpolation
        dim
        globalConnec
    end

    properties (Access = private)

        LHScells
        npnod
    end
    
    properties (Access = private)
        geometry
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSintegratorFactory();
            obj = f.create(s);
        end

    end

    methods (Access = public)
        
        function q = getQuadrature(obj)
            q = obj.quadrature;
        end
 
    end

    methods (Access = protected)
     
        function init(obj,cParams)
            obj.mesh          = cParams.mesh;
            obj.globalConnec  = cParams.globalConnec;
            obj.npnod         = cParams.npnod;
            obj.dim           = cParams.dim;
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
    
    methods (Access = private)
        
        %% LHSintegrator_triangle
      
        % Element_Elastic
        function createPrincipalDirection(obj, pdim)
            s.eigenValueComputer.type = 'PRECOMPUTED';
            s.type = pdim;
            p = PrincipalDirectionComputer.create(s);
            obj.principalDirectionComputer = p;
        end
        
        function [dir,str] = computePrincipalStressDirection(obj,tensor)
            obj.principalDirectionComputer.compute(tensor);
            dir = obj.principalDirectionComputer.direction;
            str = obj.principalDirectionComputer.principalStress;
        end

    end
    
end