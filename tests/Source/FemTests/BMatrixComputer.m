classdef BMatrixComputer < handle

    properties (Access = private)
        dim
        geometry
        globalConnec
    end
    
    methods (Access = public)

        function obj = BMatrixComputer(cParams)
            obj.init(cParams);
        end

        function Btot = compute(obj)
           Bmatrix  = obj.computeB_InMatrixForm();
           Btot     = obj.computeBtot(Bmatrix);
        end

    end

    methods (Access = private)
        
        function init(obj, cParams)
            obj.dim          = cParams.dim;
            obj.geometry     = cParams.geometry;
            obj.globalConnec = cParams.globalConnec;
        end

        function B = computeBmat(obj,igaus)
           ndim = obj.dim.ndim;
           switch ndim
               case 2
                   B = obj.computeB2D(igaus);
               case 3
                   B = obj.computeB3D(igaus);
           end
       end

       function B = computeB2D(obj,igaus)
           d = obj.dim;
           nstre          = d.nstre;
           nnode          = d.nnode;
           nelem          = d.nelem;
           nunkn          = d.nunkn;
           ndofPerElement = d.ndofPerElement;
           cartd = obj.geometry.cartd;
           B = zeros(nstre,ndofPerElement,nelem);
           for i = 1:nnode
               j = nunkn*(i-1)+1;
               B(1,j,:)   = cartd(1,i,:,igaus);
               B(2,j+1,:) = cartd(2,i,:,igaus);
               B(3,j,:)   = cartd(2,i,:,igaus);
               B(3,j+1,:) = cartd(1,i,:,igaus);
           end
       end

       function [B] = computeB3D(obj,igaus)
           d     = obj.dim;
           cartd = obj.geometry.cartd;
           B = zeros(d.nstre,d.ndofPerElement,d.nelem);
           for inode=1:d.nnode
               j = d.nunkn*(inode-1)+1;
               % associated to normal strains
               B(1,j,:)   = cartd(1,inode,:,igaus);
               B(2,j+1,:) = cartd(2,inode,:,igaus);
               B(3,j+2,:) = cartd(3,inode,:,igaus);
               % associated to shear strain, gamma12
               B(4,j,:)   = cartd(2,inode,:,igaus);
               B(4,j+1,:) = cartd(1,inode,:,igaus);
               % associated to shear strain, gamma13
               B(5,j,:)   = cartd(3,inode,:,igaus);
               B(5,j+2,:) = cartd(1,inode,:,igaus);
               % associated to shear strain, gamma23
               B(6,j+1,:) = cartd(3,inode,:,igaus);
               B(6,j+2,:) = cartd(2,inode,:,igaus);
           end
       end

       function Bmatrix = computeB_InMatrixForm(obj)
           d = obj.dim;
           Bmatrix = zeros(d.ngaus*d.nelem*d.nstre,d.ndofPerElement);
           for igaus = 1:d.ngaus
               unitaryIndex = false(d.ngaus*d.nstre,1);
               pos = d.nstre*(igaus-1) + 1 : d.nstre*(igaus);
               unitaryIndex(pos) = true;
               index = repmat(unitaryIndex,d.nelem,1);
               Bmat = obj.computeBmat(igaus);
               Bpermuted = permute(Bmat,[1 3 2]);
               Bshif = reshape(Bpermuted,d.nelem*d.nstre,d.ndofPerElement);
               Bmatrix(index,:) = Bshif;
           end
       end
       
       function Bt = computeBtot(obj, Bfull)
           d = obj.dim;
           ngaus = d.ngaus;
           nstre = d.nstre;
           nunkn = d.nunkn;
           ntot  = d.nt;
           ndofGlob = max(max(obj.globalConnec))*nunkn;
           Bt = sparse(ntot,ndofGlob);
           for idof=1:d.ndofPerElement
               GlobalDofs = obj.transformLocal2Global(idof);
               dofs = repmat(GlobalDofs',ngaus*nstre,1);
               dofs = dofs(:);
               Bt = Bt + sparse(1:ntot,dofs,Bfull(:,idof),ntot,ndofGlob);
           end
       end

       function GlobalDofs = transformLocal2Global(obj,LocalDof)
           d = obj.dim;
           connec = obj.globalConnec;
           nodesInElement   = reshape(repmat(1:d.nnode,d.nunkn,1),1,[]);
           vectorDimensions = repmat(1:d.nunkn,1,d.nnode);
           localNode        = nodesInElement(LocalDof);
           vectorDimension  = vectorDimensions(LocalDof);
           GlobalDofs       = d.nunkn*(connec(:,localNode)-1) + vectorDimension;
       end

    end
end

