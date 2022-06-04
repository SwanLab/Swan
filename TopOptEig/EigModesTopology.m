classdef EigModesTopology < handle
    
    properties (Access = public)
        V
        v1
        v2
        D
        mode1
        mode2
        xReg
    end
    
    properties (Access = private)
        designVariable
        mesh
        dim
        material
        stiffnessMatrixComputer
        massMatrixComputer
        dofConnec
        Msmooth
        filter
        iter
        elemStiff
        elemMass
        homogenizedVarComputer
        matInterpolation
        C
        dC
        Cref
        dCref
    end
    
    properties (Access = private)
        freeDOFs
        lambda
    end
    
    methods (Access = public)
        
        function obj = EigModesTopology(cParams)
            obj.init(cParams);
            obj.createBoundaryConditions();
            obj.createLHS();
            obj.iter = 0;
        end
        
        function fx = provideFunction(obj, homogenizedVariableComputer)
            obj.homogenizedVarComputer = homogenizedVariableComputer;
            obj.computeEigenModesAndValues();
            obj.lambda = obj.computeLambda();
            gamma = obj.designVariable.getFirstEigenMode();
            fx = gamma-obj.lambda(1);
        end
        
        function grad = provideDerivative(obj)
            obj.iter = obj.iter + 1;
            obj.reorderModes(obj.lambda,obj.V,obj.D);
            if obj.iter == 1
            obj.elemStiff = obj.stiffnessMatrixComputer.elemStiff;
            obj.elemMass = obj.massMatrixComputer.elemMass;
            end
            eigV1 = obj.D(1,1);
            eigV2 = obj.D(2,2);
            difEigs = abs(eigV2-eigV1);
            %if difEigs > 1
                dfdx = obj.computeSimpleEig(eigV1);
            %else
                %dfdx = obj.computeDoubleEig(Belem,x);
            %end
            % dfdx(2,nElem+1) = 1;
            g = obj.filterGradient(dfdx);
            grad = g; % (eigNum,:);
            grad(end+1,1) = 1;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.mesh           = cParams.mesh;
            obj.dim            = cParams.dim;
            obj.filter         = cParams.filter;
            obj.Msmooth        = cParams.Msmooth;
        end

        function createBoundaryConditions(obj)
            d = obj.dim;
            FixNod = obj.computeFixedNodes();
            FixDof = obj.computeFixedDOFs(FixNod);
            dofs = 1:d.ndofs;
            free  = setdiff(dofs,FixDof);
            obj.freeDOFs = free;
        end

        function FN = computeFixedNodes(obj)
            coorX = obj.mesh.coord(:,1);
            tol = 0.001;
            minX = min(coorX) + tol;
            maxX = max(coorX) - tol;
            FN1 = find(coorX < minX);
            FN2 = find(coorX > maxX);
            FN = [FN1 ;FN2];
        end

        function FixDof = computeFixedDOFs(obj, FixNod)
            d = obj.dim;
            nDOFn = d.ndofs/d.nnodes;
            FD = zeros(nDOFn*length(FixNod),1);
            for i = 1: length(FixNod)
               FD(2*i-1,1) = FixNod(i)*2-1;
               FD(2*i,1) = FixNod(i)*2;
            end
            FixDof = FD';
        end

        function createLHS(obj)
            mI = obj.createMaterialInterpolation();
            obj.createMaterial(mI);
            int = Interpolation.create(obj.designVariable.mesh,'LINEAR');
            quad = Quadrature.set(obj.designVariable.mesh.type);
            quad.computeQuadrature('QUADRATIC');
            sM.quadrature = quad;
            sM.quadType     = 'QUADRATIC';
            int.computeShapeDeriv(quad.posgp);
            sM.interpolation  = obj.mesh.interpolation;
            sS.type = 'ElasticStiffnessMatrix';
            sS.dim = obj.dim;
            sS.mesh = obj.mesh;
            sS.globalConnec   = obj.mesh.connec;
            sS.freeNodes      = obj.freeDOFs;
            sS.material       = obj.material;
            sS.interpolation  = obj.mesh.interpolation;

            obj.stiffnessMatrixComputer = LHSintegrator.create(sS);

            int = Interpolation.create(obj.mesh,'LINEAR');
            quad = Quadrature.set(obj.mesh.type);
            quad.computeQuadrature('QUADRATICMASS');
            sM.quadrature = quad;
            sM.quadType     = 'QUADRATICMASS';
            int.computeShapeDeriv(quad.posgp);
            sM.interpolation  = int;
            sM.type         = 'MassMatrix';
            sM.dim          = obj.dim;
            sM.mesh         = obj.mesh;
            sM.globalConnec = obj.mesh.connec;
            sM.freeNodes    = obj.freeDOFs;
            obj.massMatrixComputer = LHSintegrator.create(sM);

            obj.createDOFsConnec();
        end

        function createDOFsConnec(obj)
            connec  = obj.mesh.connec;
            ndimf   = obj.dim.ndimf;
            nnodeEl = size(connec, 2); 
            ndofsEl = nnodeEl * ndimf; 
            dofsElem  = zeros(ndofsEl,size(connec,1));
            for inode = 1:nnodeEl
                for iunkn = 1:ndimf
                    idofElem   = obj.nod2dof(inode,iunkn);
                    globalNode = connec(:,inode);
                    idofGlobal = obj.nod2dof(globalNode,iunkn);
                    dofsElem(idofElem,:) = idofGlobal;
                end
            end
            obj.dofConnec = dofsElem;
        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = obj.dim.ndimf;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

        function mI = createMaterialInterpolation(obj)
            designVar = obj.designVariable;
            s.dim = '2D';
            s.typeOfMaterial = 'ISOTROPIC';
            s.interpolation  = 'SIMPALL';
            s.constitutiveProperties.rho_plus = 1;
            s.constitutiveProperties.rho_minus = 1e-3;
            s.constitutiveProperties.E_plus = 1;
            s.constitutiveProperties.E_minus = 1e-3;
            s.constitutiveProperties.nu_minus = 1/3;
            s.constitutiveProperties.nu_plus = 1/3;
            s.nElem = obj.mesh.nelem;
            m = MaterialInterpolation.create(s);
            obj.matInterpolation = m;
            rho = designVar.computeVolumeFraction();
            mI = m.computeMatProp(rho);
        end

        function createMaterial(obj,mI)
            s.ptype = 'ELASTIC';
            s.pdim  = '2D';
            s.nelem = obj.mesh.nelem;
            s.mesh  = obj.mesh;
            s.kappa = mI.kappa;
            s.mu    = mI.mu;
            mat = Material.create(s);
            mat.compute(s);
            obj.material = mat;
        end

        function dfdx = computeSimpleEig(obj,eigV1)
            d = obj.dim;
            free = obj.freeDOFs;
            % ndofe = d.ndofsElem;
            % ndofn = d.ndofsElem/d.nnodeElem;
            W = zeros(d.ndofs,1);
            W(free,1) = obj.v1;
            % W(free,2) = obj.v2;
            nElem = obj.mesh.nelem;
            DOFconnec = obj.dofConnec';
%            den = obj.densityDOF(x);
            xR = obj.xReg{1};
            for i = 1:nElem
                nodDofs = DOFconnec(i,:);
                % index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
                dxS = 3*(xR(i,1))^2;
                % dxM = 1;
                % dfdx(1,i) = -(W(nodDofs,1)'*(dxS*elemStiff(:,:,i)-eigV1*dxM*massStiff(:,:,i))*W(nodDofs,1));
                dfdx(1,i) = -(W(nodDofs,1)'*(dxS*obj.elemStiff(:,:,i)-eigV1*dxS*obj.elemMass(:,:,i))*W(nodDofs,1)); %
                % dfdx(2,i) = -dx*(W(index,2)'*Belem(:,:,i)*W(index,2));
            end
        end

        function g = filterGradient(obj,dfdx)
            g = dfdx';
            gf = zeros(size(obj.Msmooth,1),obj.designVariable.nVariables);
            for ivar = 1:obj.designVariable.nVariables
                gs = g(:,:,ivar);
                gf(:,ivar) = obj.filter.getP1fromP0(gs);
            end
            gf = obj.Msmooth*gf;
            g = gf(:);
        end

        function dfdx = computeDoubleEig(obj,Belem,x)
            d    = obj.dim;
            free = obj.freeDOFs;
            ndofe = d.ndofPerElement;
            ndofn = d.ndimField;
            nElem = obj.mesh.nelem;
            W1    = zeros(d.ndof,1);
            W2    = zeros(d.ndof,1);
            dW1   = zeros(nElem,1);
            dW2   = zeros(nElem,1);
            dW1W2 = zeros(nElem,1);
            W1(free,1) = obj.v1;
            W2(free,1) = obj.v2;
            for i=1:nElem
                index = ndofn*(i-1)+1: ndofn*(i-1)+ndofe;
                dx = 2*x(i,1);
                dW1(i,1)= dx*(W1(index,1)'*Belem(:,:,i)*W1(index,1));
                dW2(i,1)= dx*(W2(index,1)'*Belem(:,:,i)*W2(index,1));
                dW1W2(i,1)= (2*x(i,1))*(W1(index,1)'*Belem(:,:,i)*W2(index,1));
                A = [dW1(i,1) dW1W2(i,1); dW1W2(i,1) dW2(i,1)];
                [U,R] = eigs(A,2,'SM');
                S = sort(diag(R));
                dfdx(1,i) = -S(1);
                dfdx(2,i) = -S(2);
            end
        end

        function computeEigenModesAndValues(obj)
            obj.filterDesignVariable()
            obj.computeTensor();
            obj.setC(obj.C)

            K = obj.stiffnessMatrixComputer.compute();
            M = obj.massMatrixComputer.compute();
            Kfree = obj.provideFreeMatrix(K);
            Mfree = obj.provideFreeMatrix(M);
            [v,d] = eigs(Kfree,Mfree,2,'SM');
            obj.V  = v;
            obj.D  = d;
            obj.printInGID();
            % obj.computeBucklingModes();
        end

        function filterDesignVariable(obj)
            x = obj.designVariable.getDensity();
            nx = length(x)/obj.designVariable.nVariables;
            xf = cell(obj.designVariable.nVariables,1);
            for ivar = 1:obj.designVariable.nVariables
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xs = x(i0:iF);
                xf{ivar,1} = obj.filter.getP0fromP1(xs);
            end
            obj.xReg = xf;
            %xf{ivar+1} = obj.designVariable.alpha;
        end

        function computeTensor(obj)
            rho = obj.xReg{1};
            nvar = 1;
            nGaus  = size(rho,2);
            nStres = obj.material.nstre;
            obj.C  = zeros(nStres,nStres,obj.mesh.nelem,nGaus);
            obj.dC = zeros(nStres,nStres,nvar,obj.mesh.nelem,nGaus);
            for igaus = 1:nGaus
                rhoV = rho(:,igaus);
                p = obj.matInterpolation.computeMatProp(rhoV);
                obj.C(:,:,:,igaus)  = obj.computeC(p.mu,p.kappa);
                for ivar = 1:nvar
                    obj.dC(:,:,1,:,igaus) = obj.computeC(p.dmu,p.dkappa);
                end
            end
            obj.Cref = obj.C;
            obj.dCref = obj.dC;
        end

        function C = computeC(obj,mu,kappa)
            s.kappa = kappa;
            s.mu    = mu;
            obj.material.compute(s);
            C  = obj.material.C;
        end

        function setC(obj,C)
            obj.material.C = C;
        end

        function printInGID(obj)
            fileName = 'Eig';
            m = obj.mesh;
            quad = Quadrature.set(m.type);
            quad.computeQuadrature('LINEAR');
            dI = obj.createPostProcessDataBase(m,fileName);
            dI.ndim   = 2;
            dI.pdim   = '2D';
            dI.ptype  = 'ELASTIC';
            dI.name = '';
            dI.printers = 'DensityGauss'; % 'HomogenizedTensor'
            p = Postprocess('VectorField',dI);
            MODE = zeros(obj.dim.ndofs,2);
            MODE(obj.freeDOFs,1) =obj.V(:,1);
            MODE(obj.freeDOFs,2) =obj.V(:,2);
            dP.fields.u = MODE;
            dP.quad   = quad;
            iter = 0;
            p.print(iter,dP);
        end

        function d = createPostProcessDataBase(obj,mesh,fileName)
            dI.mesh    = mesh;
            dI.outFileName = fileName;
            ps = PostProcessDataBaseCreator(dI);
            d = ps.create();
        end

        function MatrixFree = provideFreeMatrix(obj,Matrix)
            free = obj.freeDOFs;
            MatrixFree = Matrix(free,free);
        end

        function reorderModes(obj,lambda,V,D)
            if lambda(1)==D(1,1)
                V1=V(:,1);
                V2=V(:,2);
            else
                V1=V(:,2);
                V2=V(:,1);
            end
            obj.v1 = V1;
            obj.v2 = V2;
        end

%         function [Mode1,Mode2] = computeBucklingModes(obj)
%             nnode = obj.mesh.nnodes;
%             Mode1=zeros(nnode,2);
%             Mode2=zeros(nnode,2);
%             v1 = obj.V(:,1);
%             v2 = obj.V(:,2);
%             for i=1 : nnode-40
%             Mode1(i+40,:) = [v1(2*i-1) v1(2*i)];
%             Mode2(i+40,:) = [v2(2*i-1) v2(2*i)];
%             end
%             obj.mode1 = Mode1; 
%             obj.mode2 = Mode2; 
%         end 
        
        function l = computeLambda(obj)
            l = sort(diag(obj.D));
        end
    end
    
end