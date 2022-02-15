classdef LHSintegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        interpolation
        dim
        globalConnec
        material
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
            obj.dim          = cParams.dim;
            obj.mesh         = cParams.mesh;
            obj.npnod        = cParams.npnod;
            obj.globalConnec = cParams.globalConnec;
            if isfield(cParams, 'material') % Compatibility with MassMatrix
                obj.material   = cParams.material;
            end
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
%             ndofs  = obj.npnod; % should be obj.dim.ndof
            ndofs  = obj.dim.ndof; % should be obj.dim.ndof
            Ae     = aElem;
            nunkn1 = obj.dim.nunkn;
            nunkn2 = obj.dim.nunkn;
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

%         function A = assembleMatrix(obj,aElem)
%             connec = obj.globalConnec;
%             nunkn     = obj.dim.nunkn;
%             nnode     = obj.dim.nnode;
%             ndofs     = obj.dim.ndof; % should be obj.dim.ndof
%             ndofsElem = obj.dim.ndofPerElement;
%             Ae     = aElem;
%             A = sparse(ndofs,ndofs);
%             for iunkn = 1:nunkn
%                 for i = 1:nnode
%                     nodeI = connec(:,i);
%                     NODEI = nunkn*(nodeI-1)+iunkn;
%                     for j = 1:nnode
%                         nodeJ = connec(:,j); %DOF = ni*(nod-1)+dir;
%                         NODEJ = nunkn*(nodeJ-1)+iunkn;
%                         a = squeeze(Ae(i,j,:));
%                         A = A + sparse(NODEI,NODEJ,a,ndofs,ndofs);
%                     end
%                 end
%             end
%         end

%         function A = assembleMatrix(obj,aElem)
%             connec = obj.globalConnec;
%             nunkn     = obj.dim.nunkn;
%             nnode     = obj.dim.nnode;
%             ndofs     = obj.dim.ndof; % should be obj.dim.ndof
%             ndofsElem = obj.dim.ndofPerElement;
%             Ae     = aElem;
%             A = sparse(ndofs,ndofs);
%             for iunkn = 1:nunkn
%                 for i = 1:nnode
%                     nodeI = connec(:,i);
%                     NODEI = nunkn*(nodeI-1)+iunkn;
%                     for j = 1:nnode
%                         nodeJ = connec(:,j); %DOF = ni*(nod-1)+dir;
%                         NODEJ = nunkn*(nodeJ-1)+iunkn;
%                         a = squeeze(Ae(i,j,:));
%                         A = A + sparse(NODEI,NODEJ,a,ndofs,ndofs);
%                     end
%                 end
%             end
%         end

    end
    
    methods (Access = private)
        
        function Aglo = nod2dof(obj,AnodGLO,ndim)
        end

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