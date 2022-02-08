classdef LHSintegrator_triangle < LHSintegrator
    
    properties(Access = private)
%         K_generator
%         Bmatrix
        geometry
    end

    methods (Access = public)

        function Kgen = computeKgenerator(obj)
            % Previously computeTriangleLHS
            % Should really become a LHS computer and avoid using K
            % generators, using shape functions instead.
           obj.createGeometry();
           Bmat = obj.computeBmat();
           connect = obj.mesh.connec;
           dvolum = obj.geometry.dvolu;
%            obj.K_generator = StiffnessMatrixGenerator(connect,Bmat,dvolum,obj.dim);
           Bmatrix = obj.computeB_InMatrixForm();
           Kgen = KGeneratorWithfullStoredB(obj.dim,connect,Bmatrix,dvolum);
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