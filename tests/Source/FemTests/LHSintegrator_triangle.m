classdef LHSintegrator_triangle < LHSintegrator
    
    properties(Access = private)
        K_generator
        Bmatrix
        geometry        
    end

    methods (Access = public)

        function LHS = computeTriangleLHS(obj)
           obj.createGeometry();
           Bmat = obj.computeBmat();
%            ngaus   = obj.quadrature.ngaus;
%            obj.dim   = obj.computeDim(ngaus);
           connect = obj.mesh.connec;
           dvolum = obj.geometry.dvolu;
           obj.K_generator = StiffnessMatrixGenerator(connect,Bmat,dvolum,obj.dim);
           obj.Bmatrix = obj.computeB_InMatrixForm();
           LHS = KGeneratorWithfullStoredB(obj.dim,connect,obj.Bmatrix,dvolum);
        end
    end

   methods (Access = private)

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
      
       function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interpolation);
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
            nnode = obj.mesh.nnode;
            nelem = obj.mesh.nelem;
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

        function Bmatrix = computeB_InMatrixForm(obj)
            ndofPerElement = obj.dim.ndofPerElement;
            Bfull = zeros(obj.quadrature.ngaus,obj.dim.nstre,ndofPerElement,obj.dim.nelem);
            Bmatrix = zeros(obj.quadrature.ngaus*obj.dim.nelem*obj.dim.nstre,ndofPerElement);
            for igaus = 1:obj.quadrature.ngaus
                unitaryIndex = false(obj.quadrature.ngaus*obj.dim.nstre,1);
                pos = obj.dim.nstre*(igaus-1) + 1 : obj.dim.nstre*(igaus) ;
                unitaryIndex(pos) = true;
                
                Index = repmat(unitaryIndex,obj.dim.nelem,1);
                Bfull(igaus,:,:,:) = obj.computeB(igaus);
                Bshif = reshape(permute(obj.computeB(igaus),[1 3 2]),obj.dim.nelem*obj.dim.nstre,ndofPerElement);
                Bmatrix(Index,:) = Bshif;
            end
        end
   end

end