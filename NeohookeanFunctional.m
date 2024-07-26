classdef NeohookeanFunctional < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        mesh
        lambda
        mu
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = NeohookeanFunctional(cParams)
            obj.init(cParams)
%             obj.mu = 10 / (2*(1 + 0.3));
%             obj.lambda = (10*0.3) / ((1 + 0.3) * (1 - 2*0.3));
        end

        function val = compute(obj, uFun)
            quad = Quadrature.create(obj.mesh, 2);
            xG = quad.posgp;
            nDimf = uFun.ndimf;
            
            [F,~] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);
            
            C = pagemtimes(Ft,F);
            trC = obj.computeTrace(C);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);

            % val = obj.mu/2*(trC - nDimf) - obj.mu*log(jac) + obj.lambda/2*(jac-1).^2;
            val = obj.mu/2*(trC - nDimf) - obj.mu*log(jac) + obj.lambda/2*(log(jac)).^2; % stanford
            val = pagemtimes(val,dV);
            val = sum(val, 'all');
        end
        

        function Fint = computeGradient(obj, uFun)
            nDimf = uFun.ndimf;
            test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            quad = Quadrature.create(obj.mesh,3);

            xG = quad.posgp;
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);
            dNdxTest  = test.evaluateCartesianDerivatives(xG);

            nNode = size(dNdxTest,2);
            nGaus = size(dNdxTest,3);
            nElem = size(dNdxTest,4);
            nDof = nDimf*nNode;
            nDim = obj.mesh.ndim;

            piola = obj.computeFirstPiola(uFun,xG);
            dofToDim = repmat(1:nDimf,[1,nNode]);
            dofToNode = repmat(1:nNode, nDimf, 1);
            dofToNode = dofToNode(:);

            fint = zeros(nDof,1,nGaus,nElem);
            fint2 = zeros(nDof,1,nGaus,nElem);


            for iNode = 1:nNode
                for iDim = 1:nDim
                    for iField = 1:nDimf
                        dNdx_ij = dNdxTest(iDim, iNode, :, :);
                        iDof = iField + nDimf*(iNode-1);
                        Pik = piola(iField,iDim,:,:);
                        fint2(iDof,1,:,:) = fint2(iDof,1,:,:) + Pik.*dNdx_ij;
                    end
                end
            end

            % GradDeltaV is not always compatible (see BCs), but we dont
            % worry about it since we reduce the matrix later on    
            for iDof = 1:nDof
                iNode = dofToNode(iDof);
                iDim  = dofToDim(iDof);
                deltav = zeros(nNode,nDimf, nGaus, nElem);
                deltav(iNode,iDim,:,:) = 1;
                GradDeltaV = pagemtimes(dNdxTest,deltav);
                fint(iDof, :,:,:) = squeeze(bsxfun(@(A,B) sum(A.*B, [1 2]), pagetranspose(piola),GradDeltaV));
            end
            err = norm(fint(:)-fint2(:))/norm(fint(:))
            fint = fint2;
            fint = fint.*dV;
            fint = squeeze(sum(fint,3));
            Fint = obj.assembleIntegrand(fint,test);
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = pagetranspose(rhsElem);
            connec = test.getConnec();
            nDofs = max(max(connec));
            nDofElem  = size(connec,2);
            f = zeros(nDofs,1);
            for idof = 1:nDofElem
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end
        end

        function hess = computeHessian(obj, uFun)
            % This is the LINEALIZED hessian (Holzapfel, 401)
            % See  Holzapfel, 396
            nDimf = uFun.ndimf;
            trial = uFun;
            test  = LagrangianFunction.create(obj.mesh, nDimf, 'P1');
            quad = Quadrature.create(obj.mesh,3);

            xG = quad.posgp;
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);

            [F,~] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F,[2 1 3 4]);
            Aneo = obj.computeTangentConstitutive(uFun,xG);
            dNdxTest  = test.evaluateCartesianDerivatives(xG);
            dNdxTrial = trial.evaluateCartesianDerivatives(xG);
            nNode = size(dNdxTrial,2);
            nGaus = size(dNdxTrial,3);
            nElem = size(dNdxTrial,4);
            nDof = nDimf*nNode;
            
            dofToDim = repmat(1:nDimf,[1,nNode]);
            dofToNode = repmat(1:nNode, nDimf, 1);
            dofToNode = dofToNode(:);

            K = zeros(nDof,nDof,nGaus,nElem);
            for iDof = 1:nDof % test dof
                iNode = dofToNode(iDof);
                iDim  = dofToDim(iDof);
                GradDeltaV = zeros(nDimf,nDimf, nGaus, nElem);
                GradDeltaV(:,iDim,:,:) = dNdxTest(:,iNode,:,:);
%                 GradDeltaV = pagemtimes(F,GradDeltaV);
                res = zeros(nDimf,nDimf,nGaus,nElem);
                for a = 1:nDimf
                    for b = 1:nDimf
                        C = squeezeParticular(Aneo(:,:,a,b,:,:),[3 4]);
                        res(a,b,:,:) = bsxfun(@(A,B) sum(A.*B, [1 2]), GradDeltaV,pagetranspose(C));
                    end
                end

                for jDof = 1:nDof % trial dof
                    jNode = dofToNode(jDof);
                    jDim  = dofToDim(jDof);

                    GradDeltaU = zeros(nDimf,nDimf, nGaus, nElem);
                    GradDeltaU(:,jDim,:,:) = dNdxTrial(:,jNode,:,:);
%                     GradDeltaU = pagemtimes(GradDeltaU,Ft);

%                     Kgeo = bsxfun(@(A,B) sum(A.*B, [1 2]), piola, GradDeltaU);

                   Ktan = bsxfun(@(A,B) sum(A.*B, [1 2]), pagetranspose(res), GradDeltaU);
                    K(iDof,jDof,:,:) = Ktan;
                end
            end
            K = K.*dV;
            K = squeeze(sum(K,3));

            % Assembly
            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            hess = assembler.assemble(K, test, trial);
        end


    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.mesh   = cParams.mesh;
        end

        function piola = computeFirstPiola(obj,uFun,xG)
            [F,~] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);

            nGaus = size(xG,2);
            nElem = size(Ft,4);

            invFt = MatrixVectorizedInverter.computeInverse(Ft);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            jac = reshape(jac, [1 1 nGaus nElem]);
            piola = obj.mu*(F-invFt) + obj.lambda*log(jac).*invFt;
%             piola = obj.mu*(F-invFt) + obj.lambda/2*(jac.^2-1).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun,xG)
            [F,I33] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);
            
            nGaus = size(xG,2);
            nElem = size(Ft,4);

            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = MatrixVectorizedInverter.computeInverse(Ft);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            jac = reshape(jac, [1 1 1 1 nGaus nElem]);
            logJac = log(jac);

            Aneo = obj.lambda*obj.outerProduct(invFt, invFt) + ...
                obj.mu*obj.kron_topF(I33,I33) + ...
                (obj.mu-obj.lambda*logJac).*obj.kron_botF(invFt, invF);
%             quad = Quadrature.create(obj.mesh,3);
%             xG = quad.posgp;

%             [F,~] = obj.computeDeformationGradient(uFun, xG);
%             Ft = permute(F, [2 1 3 4]);
%             nF    = size(F,1);
%             nGaus = size(F,3);
%             nElem = size(F,4);
%             A = zeros(nF,nF,nF,nF,nGaus,nElem);
%                 for i = 1:nF
%                 for j = 1:nF
%                 for k = 1:nF
%                 for l = 1:nF
%                     first_term(1,1,1,1,:,:) = squeeze(obj.lambda*invF(j,i,:,:).*invF(l,k,:,:));
%                     secnd_term(1,1,1,1,:,:) = squeeze(obj.mu*I33(i,k,:,:).*I33(j,l,:,:));
%                     trd(1,1,:,:,:,:) = invF(l,i,:,:).*invF(j,k,:,:);
%                     third_term(1,1,1,1,:,:) = squeeze((obj.mu - obj.lambda*log(jac)).*trd);
%                     A(i,j,k,l,:,:) = first_term + secnd_term + third_term;
% %                     A(i,j,k,l,:,:) = first_term;
%                 end
%                 end
%                 end
%                 end
%             norm(Aneo(:)-A(:))
%             Aneo = A;
        end

        function [F,I33] = computeDeformationGradient(obj, uFun, xG)
            nPoints  = size(xG,2);
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;

            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
            GradU = permute(GradU, [2 1 3 4]);
%             GradU = [0.0 0.0 0.0; -3.415063509461096 -0.24999999999999956 -0.4330127018922192; 0.9150635094610968 0.43301270189221924 -0.24999999999999994];

            I33 = obj.createIdentityMatrix(size(GradU));

            F = I33 + GradU;
        end

        function trC = computeTrace(obj,C)
            trC = zeros(1,1,size(C,3),size(C,4));
            for i = 1:size(C,1)
                trC = trC + C(i,i,:,:,:);
            end
        end

        function I = createIdentityMatrix(obj,sz)
            I = zeros(sz);
            for i = 1:sz(1)
                I(i,i,:,:) = 1;
            end
        end

        % Operators
        function C = double_dot(obj,A, B)
            assert(~isvector(A) && ~isvector(B))
            idx = max(0, ndims(A) - 1);
            B_t = permute(B, circshift(1:ndims(A) + ndims(B), [0, idx - 1]));
            C = squeeze(sum(squeeze(sum(bsxfun(@times, A, B_t), idx)), idx));
        end

        function C = outerProduct(obj, A, B)  % version 1
            C = zeros([size(A,1), size(A,2),size(B)]);
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(B,1)
                        for l = 1:size(B,2)
                            C(i,j,k,l,:,:) = A(i,j,:,:).*B(k,l,:,:);
                        end
                    end
                end
            end
        end

%         function C= kron_topF(obj, A,B)
%             C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
%             for i = 1:size(A,1)
%                 for k = 1:size(B,1)
%                     C(i,:,k,:,:,:) = pagemtimes( A(i,k,:,:), B);
%                 end
%             end
%         end

        function C= kron_topF(obj,A,B)
    %             C = obj.kron_topF(A,B);
    %             C = pagetranspose(C);
            C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
            for i = 1:size(A,1)
                for k = 1:size(A,2)
                    for j = 1:size(B,1)
                        for l = 1:size(B,2)
                            C(i,j,k,l,:,:) = A(i,k,:,:).*B(j,l,:,:);
                        end
                    end
                end
            end
        end

        function C= kron_botF(obj,A,B)
    %             C = obj.kron_topF(A,B);
    %             C = pagetranspose(C);
            C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
            for i = 1:size(A,1)
                for j = 1:size(A,2)
                    for k = 1:size(B,1)
                        for l = 1:size(B,2)
                            C(i,j,k,l,:,:) = A(i,l,:,:).*B(j,k,:,:);
                        end
                    end
                end
            end
        end
        
    end
    
end