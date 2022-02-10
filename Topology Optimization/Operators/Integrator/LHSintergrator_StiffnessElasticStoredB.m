classdef LHSintergrator_StiffnessElasticStoredB < LHSintegrator
  
    properties (Access = public)
        K
    end

    properties (Access = private)
%         dofsPerElement
%         nodesInElement
%         VectorDimensions
%         Bfull
%         nunkn
%         dvolum
        ndofGlobal
        nt
    end

    properties (Access = private)
        material
        geometry
        Btot
        nodesInElement
        VectorDimensions
    end

    methods (Access = public)
        
        function obj = LHSintergrator_StiffnessElasticStoredB(cParams)
            obj.init(cParams)
            obj.material = cParams.material;
            obj.createQuadrature();
            obj.createInterpolation();

            % B calculations from Element
            obj.createGeometry();
            Bmat = obj.computeBmat();
            Bmatrix = obj.computeB_InMatrixForm();

            % B calculations from KGenerator
            d = obj.dim;
            obj.nt = d.ngaus*d.nelem*d.nstre;
            obj.ndofGlobal = max(max(obj.globalConnec))*d.nunkn;
            obj.nodesInElement   = reshape(repmat(1:d.nnode,d.nunkn,1),1,[]);
            obj.VectorDimensions = repmat(1:d.nunkn,1,d.nnode);
            obj.Btot = obj.computeBtot(Bmatrix);
        end

        function LHS = compute(obj)
%             lhs = obj.computeElementalLHS();
%             LHS = obj.assembleMatrix(lhs);
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

       function createGeometry(obj)
           s.mesh = obj.mesh;
           obj.geometry = Geometry.create(s);
           obj.geometry.computeGeometry(obj.quadrature,obj.interpolation);
       end

       function Bmatrix = computeB_InMatrixForm(obj)
           d = obj.dim;
           Bfull = zeros(d.ngaus,d.nstre,d.ndofPerElement,d.nelem);
           Bmatrix = zeros(d.ngaus*d.nelem*d.nstre,d.ndofPerElement);
           for igaus = 1:d.ngaus
               unitaryIndex = false(d.ngaus*d.nstre,1);
               pos = d.nstre*(igaus-1) + 1 : d.nstre*(igaus) ;
               unitaryIndex(pos) = true;
               
               Index = repmat(unitaryIndex,d.nelem,1);
               Bfull(igaus,:,:,:) = obj.computeB(igaus);
               Bshif = reshape(permute(obj.computeB(igaus),[1 3 2]),d.nelem*d.nstre,d.ndofPerElement);
               Bmatrix(Index,:) = Bshif;
           end
       end

       function Bmat = computeBmat(obj)
            ngaus = obj.quadrature.ngaus;
            nelem = obj.mesh.nelem;
            nstre = obj.dim.nstre;
            ndofPerElement = obj.dim.ndofPerElement;
            Bmat = zeros(ngaus,nstre,ndofPerElement,nelem);
            for igaus = 1:ngaus
                Bmat(igaus,:,:,:) = obj.computeB(igaus);
            end
       end

       function B = computeB(obj,igaus)
           ndim = obj.dim.ndim;
           switch ndim
               case 2
                   B = obj.computeB2D(igaus);
               case 3
                   B = obj.computeB3D(igaus);
           end
       end

       function B = computeB2D(obj,igaus)
            nstre = obj.dim.nstre;
            nnode = obj.dim.nnode;
            nelem = obj.dim.nelem;
            nunkn = obj.dim.nunkn; 
            ndofPerElement = obj.dim.ndofPerElement;
            B = zeros(nstre,ndofPerElement,nelem);
            for i = 1:nnode
                j = nunkn*(i-1)+1;
                B(1,j,:)  = obj.geometry.cartd(1,i,:,igaus);
                B(2,j+1,:)= obj.geometry.cartd(2,i,:,igaus);
                B(3,j,:)  = obj.geometry.cartd(2,i,:,igaus);
                B(3,j+1,:)= obj.geometry.cartd(1,i,:,igaus);
            end
       end

       function [B] = computeB3D(obj,igaus)
           d = obj.dim;
            B = zeros(d.nstre,d.ndofPerElement,d.nelem);
            for inode=1:d.nnode
                j = d.nunkn*(inode-1)+1;
                % associated to normal strains
                B(1,j,:) = obj.geometry.cartd(1,inode,:,igaus);
                B(2,j+1,:) = obj.geometry.cartd(2,inode,:,igaus);
                B(3,j+2,:) = obj.geometry.cartd(3,inode,:,igaus);
                % associated to shear strain, gamma12
                B(4,j,:) = obj.geometry.cartd(2,inode,:,igaus);
                B(4,j+1,:) = obj.geometry.cartd(1,inode,:,igaus);
                % associated to shear strain, gamma13
                B(5,j,:) = obj.geometry.cartd(3,inode,:,igaus);
                B(5,j+2,:) = obj.geometry.cartd(1,inode,:,igaus);
                % associated to shear strain, gamma23
                B(6,j+1,:) = obj.geometry.cartd(3,inode,:,igaus);
                B(6,j+2,:) = obj.geometry.cartd(2,inode,:,igaus);
            end
       end

       function CmatTot = computeCmatBlockDiagonal(obj)
           Cmat = obj.material.C;
           CmatTot = sparse(obj.nt,obj.nt);
           dvol = obj.geometry.dvolu;
           for istre = 1:obj.dim.nstre
               for jstre = 1:obj.dim.nstre
                   for igaus = 1:obj.dim.ngaus
                       posI = (istre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : obj.nt ;
                       posJ = (jstre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : obj.nt ;
                       
                       Ct = squeeze(Cmat(istre,jstre,:,igaus)).*dvol(:,igaus);
                       CmatTot = CmatTot + sparse(posI,posJ,Ct,obj.nt,obj.nt);
                   end
               end
           end
       end
       
       function Bt = computeBtot(obj, Bfull)
           Bt = sparse(obj.nt,obj.ndofGlobal);
           d = obj.dim;
           for idof=1:d.ndofPerElement
               GlobalDofs = obj.transformLocal2Global(idof);
               dofs = repmat(GlobalDofs',d.ngaus*d.nstre,1);
               dofs = dofs(:);
               Bt = Bt + sparse(1:obj.nt,dofs,Bfull(:,idof),obj.nt,obj.ndofGlobal);
           end
%            obj.Btot = Bt;
       end


       function GlobalDofs = transformLocal2Global(obj,LocalDof)
           LocalNode        = obj.nodesInElement(LocalDof);
           VectorDimension  = obj.VectorDimensions(LocalDof);
           GlobalDofs       = obj.dim.nunkn*(obj.globalConnec(:,LocalNode)-1) + VectorDimension;
       end

       function K = computeStiffness(obj,CmatTot)
           B = obj.Btot;
           CB = CmatTot*B;
           K = B'*CB;
       end

   end
    
end