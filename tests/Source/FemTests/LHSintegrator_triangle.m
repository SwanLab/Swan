classdef LHSintegrator_triangle < LHSintegrator
    
    properties(Access = private)
        K_generator
        Bmatrix
        StiffnessMatrix
        geometry
        dim
        material
    end

    properties (Access = public)
        bcApplier
        Kred
        dof
    end

    methods (Access = public)
        function copiat(obj)
            obj.initcopiat();
        end
    end
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
            dShape2 = obj.interpolation.deriv;
            dShape = obj.computeGradient();
            dvolu  = obj.mesh.computeDvolume(obj.quadrature);
            ngaus  = obj.quadrature.ngaus;
            nelem  = obj.mesh.nelem;
            nnode  = obj.mesh.nnode;
            lhs = zeros(nnode,nnode,nelem);
            for igaus = 1:ngaus
                dv(1,1,:) = dvolu(igaus,:);
                for iNode = 1:nnode
                   for jNode = 1:nnode 
                      dNi = dShape(:,iNode,:,igaus);
                      dNj = dShape(:,jNode,:,igaus);
                      dNidNj = sum(dNi.*dNj,1);
                      lhs(iNode,jNode,:) = lhs(iNode,jNode,:) + dNidNj.*dv;
                   end
                end
            end
        end
        
   end
    
   methods (Access = private)

       function grad = computeGradient(obj)
            q = Quadrature.set(obj.mesh.type);
            q.computeQuadrature('LINEAR');
            m.type = obj.mesh.type;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.quadrature.posgp);
            s.mesh = obj.mesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,int);
            grad = g.cartd;
       end

       function initcopiat(obj)
           obj.geometria();
            Bmat = obj.computeBmat();
            ngaus   = obj.quadrature.ngaus;
            dimen   = obj.computeDim(ngaus);
            obj.dim = dimen;
            connect = obj.mesh.connec;
            dvolum = obj.geometry.dvolu;
            obj.K_generator = StiffnessMatrixGenerator(connect,Bmat,dvolum,dimen);
            obj.Bmatrix = obj.computeB_InMatrixForm();
            obj.StiffnessMatrix = KGeneratorWithfullStoredB(dimen,connect,obj.Bmatrix,dvolum);
            dEps = obj.computedEps();
            K = obj.computeStiffnessMatrixSYM();
            bca = obj.createBCApplier();
            obj.bcApplier = bca;
            obj.Kred = bca.fullToReducedMatrix(K);
       end

        function bcApplier = createBCApplier(obj)
            nfields = 1;
            interp{1} = obj.interpolation; % legacy, apparently
            obj.dof = DOF_Elastic(obj.fileName,obj.mesh,obj.problemData.pdim,nfields,interp);
            cParams.nfields = nfields; % hmmmm, apparently it is so
            cParams.dof = obj.dof;
            cParams.scale = obj.problemData.scale;
            cParams.type = 'Dirichlet'; % defined in Element
            bcApplier = BoundaryConditionsApplier.create(cParams);
        end

        % copiat de Element_Elastic
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
       %copiat de ElasticDim
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

       %copiat de ElasticDim
       function dEps = computedEps(obj)
            dvolum = obj.geometry.dvolu';
            nv     = obj.dim.nnode*2;
            dEps = zeros(nv,obj.quadrature.ngaus,obj.dim.nstre,obj.dim.nelem);
            for igaus = 1:obj.quadrature.ngaus
               B  = obj.computeB(igaus); 
               Bm = permute(B,[2 1 3]);
               dvG(1,1,:) = squeeze(dvolum(igaus,:));
               dvGm = repmat(dvG,nv,obj.dim.nstre,1); % 3 hardcoded
               dEps(:,igaus,:,:) = Bm.*dvGm;
            end
       end
      
       function geometria(obj)
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
       %copiat de Element_Elastic
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

        function fesMaterial(obj)
            s.ptype = obj.problemData.ptype;
            s.pdim  = obj.problemData.pdim;
            s.nelem = obj.mesh.nelem;
            s.geometry = obj.geometry;
            s.mesh  = obj.mesh;            
            obj.material = Material.create(s);
        end

        function K = computeStiffnessMatrixSYM(obj)
            obj.fesMaterial();
            obj.computeC();
            obj.StiffnessMatrix.compute(obj.material.C);
            K = obj.StiffnessMatrix.K;
        end

        function computeC(obj) %IsotropicElasticMaterial
            I = ones(obj.mesh.nelem,obj.quadrature.ngaus);
            kappa = .9107*I;
            mu    = .3446*I;
            nElem = size(mu,1);
            nGaus = size(mu,2);
            m = mu;
            l = kappa - mu;
            C = zeros(obj.dim.nstre,obj.dim.nstre,nElem,nGaus);             
            C(1,1,:,:)= 2*m+l;
            C(1,2,:,:)= l;
            C(2,1,:,:)= l;
            C(2,2,:,:)= 2*m+l;
            C(3,3,:,:)= m;
            obj.material.C = C;
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