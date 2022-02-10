classdef LHSintergrator_StiffnessElasticStoredB < LHSintegrator

    properties (Access = private)
        nt
        material
        geometry
        Btot
    end

    methods (Access = public)
        
        function obj = LHSintergrator_StiffnessElasticStoredB(cParams)
            obj.init(cParams)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.initOwn(cParams);
            obj.createGeometry();
            obj.computeB();
        end

        function LHS = compute(obj)
            CmatTot = obj.computeCmatBlockDiagonal();
            LHS = obj.computeStiffness(CmatTot);
        end
        
        function setMaterialC(obj, Cmat)
            obj.material.C = Cmat;
        end
    end
    
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
        end

        function lhs = assembleMatrix(obj)
        end
        
   end
    
   methods (Access = private)
       
       function initOwn(obj, cParams)
            obj.material = cParams.material;
            d            = obj.dim;
            obj.nt       = d.ngaus*d.nelem*d.nstre;
       end

       function createGeometry(obj)
           s.mesh = obj.mesh;
           obj.geometry = Geometry.create(s);
           obj.geometry.computeGeometry(obj.quadrature,obj.interpolation);
       end

       function computeB(obj)
           Bmatrix  = obj.computeB_InMatrixForm();
           obj.Btot = obj.computeBtot(Bmatrix);
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
           d = obj.dim;
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
       
       function Bt = computeBtot(obj, Bfull)
           ndofGlob = max(max(obj.globalConnec))*obj.dim.nunkn;
           Bt = sparse(obj.nt,ndofGlob);
           d = obj.dim;
           ngaus = d.ngaus;
           nstre = d.nstre;
           ntot  = obj.nt;
           for idof=1:d.ndofPerElement
               GlobalDofs = obj.transformLocal2Global(idof);
               dofs = repmat(GlobalDofs',ngaus*nstre,1);
               dofs = dofs(:);
               Bt = Bt + sparse(1:ntot,dofs,Bfull(:,idof),ntot,ndofGlob);
           end
       end

       function CmatTot = computeCmatBlockDiagonal(obj)
           nstre = obj.dim.nstre;
           ngaus = obj.dim.ngaus;
           ntot  = obj.nt;
           Cmat    = obj.material.C;
           CmatTot = sparse(ntot,ntot);
           dvol = obj.geometry.dvolu;
           for istre = 1:nstre
               for jstre = 1:nstre
                   for igaus = 1:ngaus
                       posI = (istre)+(nstre)*(igaus-1) : ngaus*nstre : ntot;
                       posJ = (jstre)+(nstre)*(igaus-1) : ngaus*nstre : ntot ;
                       
                       Ct = squeeze(Cmat(istre,jstre,:,igaus)).*dvol(:,igaus);
                       CmatTot = CmatTot + sparse(posI,posJ,Ct,ntot,ntot);
                   end
               end
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

       function K = computeStiffness(obj,CmatTot)
           B = obj.Btot;
           CB = CmatTot*B;
           K = B'*CB;
       end

   end
    
end