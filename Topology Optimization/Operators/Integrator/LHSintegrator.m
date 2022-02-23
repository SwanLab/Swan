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
            connec    = obj.globalConnec;
            dofConnec = obj.computeDOFconnec();
            ndofs  = obj.dim.ndof;
            nunkn  = obj.dim.nunkn;
%             nnode  = obj.dim.nnode;
            nnode  = size(connec,2); % pending review why TopOptTests take incorrect nnode
            Ae     = aElem;
            A = sparse(ndofs,ndofs);
            for i = 1:nnode*nunkn
                dofsI = dofConnec(:,i);
                for j = 1:nnode*nunkn
                    dofsJ = dofConnec(:,j);
                    a = squeeze(Ae(i,j,:));
                    Aadd = obj.computeAaddBySparse(a, dofsI, dofsJ);
%                     Aadd = obj.computeAaddByAccumarray(a, dofsI, dofsJ);
                    A = A + Aadd;
                end
            end
        end

    end
    
    methods (Access = private)

        function Cadd = computeAaddBySparse(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndof;
           Cadd = sparse(dofsI,dofsJ,a,ndofs,ndofs);
        end

        function Cadd = computeAaddByAccumarray(obj,a, dofsI, dofsJ)
           d = obj.dim;
           ndofs = d.ndof;
           index = [dofsI', dofsJ'];
           Cadd = accumarray(index,a,[ndofs ndofs],[],[],true);
        end

        function dof_elem = computeDOFconnec(obj)
            connec = obj.globalConnec;
            nunkn  = obj.dim.nunkn;
            nnode  = size(connec,2);
            dof_elem  = zeros(nnode*nunkn,size(connec,1));
            for inode = 1:nnode
                for iunkn = 1:nunkn
                    idof_elem = obj.nod2dof(inode,iunkn);
                    global_node = connec(:,inode);
                    idof_global = obj.nod2dof(global_node,iunkn);
                    dof_elem(idof_elem,:) = idof_global;
                end
            end
            dof_elem = dof_elem';
        end

        function idof = nod2dof(obj, inode, iunkn)
            nunkn = obj.dim.nunkn;
            idof(:,1)= nunkn*(inode - 1) + iunkn;
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