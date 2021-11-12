classdef LHSintegrator_triangle < LHSintegrator
    
    properties(Access = private)
        K_generator
        Bmatrix
        geometry
        dim
        material
    end

    properties (Access = public)
        Kred
        dof
        StiffnessMatrix
    end

    methods (Access = public)

        function computeTriangleLHS(obj)
           obj.createGeometry();
           Bmat = obj.computeBmat();
           ngaus   = obj.quadrature.ngaus;
           obj.dim   = obj.computeDim(ngaus);
           connect = obj.mesh.connec;
           dvolum = obj.geometry.dvolu;
           obj.K_generator = StiffnessMatrixGenerator(connect,Bmat,dvolum,obj.dim);
           obj.Bmatrix = obj.computeB_InMatrixForm();
           obj.StiffnessMatrix = KGeneratorWithfullStoredB(obj.dim,connect,obj.Bmatrix,dvolum);
        end
    end

   methods (Access = private)

       % Element_Elastic
       function dim = computeDim(obj,ngaus)
           dim                = DimensionVariables();
           dim.nnode          = obj.mesh.nnode;
           dim.nunkn          = 2;
           dim.nstre          = 3;
           dim.ndof           = obj.npnod*dim.nunkn;
           dim.nelem          = obj.mesh.nelem;
           dim.ndofPerElement = dim.nnode*dim.nunkn;
           dim.ngaus          = ngaus;
           dim.nentries       = dim.nelem*(dim.ndofPerElement)^2;
       end
       
       % ElasticDim
       function Bmat = computeBmat(obj)
            ngaus = obj.quadrature.ngaus;
            nnode = obj.interpolation.nnode;
            nunkn = 2; %hardcoded
            nelem = obj.mesh.nelem;
            nstre = 3; %hardcoded
            Bmat = zeros(ngaus,nstre,nnode*nunkn,nelem); %hardcoded
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
            nstre = 3;
            nnode  = obj.mesh.nnode;
            nelem  = obj.mesh.nelem;
            B = zeros(nstre,nnode*2,nelem); % HARDCODED
            for i = 1:nnode
                j = 2*(i-1)+1; % HARDCODED
                B(1,j,:)  = obj.geometry.cartd(1,i,:,igaus);
                B(2,j+1,:)= obj.geometry.cartd(2,i,:,igaus);
                B(3,j,:)  = obj.geometry.cartd(2,i,:,igaus);
                B(3,j+1,:)= obj.geometry.cartd(1,i,:,igaus);
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