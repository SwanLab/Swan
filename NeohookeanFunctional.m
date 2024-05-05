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

            nPoints  = quad.ngaus;
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;
            GradU = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);

            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1;
            I33(2,2,:,:) = 1;
            I33(3,3,:,:) = 1;

            F = I33 + GradU; % deformation gradient
            F = permute(F, [2 1 3 4]);
            Ft = permute(F, [2 1 3 4]);
            
            C = pagemtimes(Ft,F);
            trC = C(1,1,:,:) +C(2,2,:,:)+C(3,3,:,:);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);

            % val = obj.mu/2*(trC - 3) - obj.mu*log(jac) + obj.lambda/2*(jac-1).^2;
            val = obj.mu/2*(trC - 3) - obj.mu*log(jac) + obj.lambda/2*(log(jac)).^2; % stanford
            val = pagemtimes(val,dV);
            val = sum(val, 'all');
        end

        function Fint = computeInternalForces(obj, uFun)
            test = LagrangianFunction.create(obj.mesh, obj.mesh.ndim, 'P1');
            quad = Quadrature.create(obj.mesh,2);

            xG = quad.posgp;
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);
            dNdxTest  = test.evaluateCartesianDerivatives(xG);

            nNode = size(dNdxTest,2);
            nGaus = size(dNdxTest,3);
            nElem = size(dNdxTest,4);

            piola = obj.computeFirstPiola(uFun,xG);
            dofToDim = repmat(1:3,[1,nNode]);
            dofToNode = repmat(1:nNode,[1,3]);

            fint = zeros(3*nNode,1,nGaus,nElem);
            for iDof = 1:24
                iNode = dofToNode(iDof);
                iDim  = dofToDim(iDof);
                GradDeltaV = zeros(3,3, nGaus, nElem);
                GradDeltaV(iDim,:,:,:) = dNdxTest(:,iNode,:,:);
                GradDeltaV = permute(GradDeltaV, [2 1 3 4]);
                fint(iDof, :,:,:) = squeeze(bsxfun(@(A,B) sum(A.*B, [1 2]), GradDeltaV,piola));
            end
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
            trial = uFun;
            test  = LagrangianFunction.create(obj.mesh, 3, 'P1');
            quad = Quadrature.create(obj.mesh,2);

            xG = quad.posgp;
            dV(1,1,:,:) = obj.mesh.computeDvolume(quad);

            Ctan = obj.computeTangentConstitutive(uFun,xG);

            dNdxTest  = test.evaluateCartesianDerivatives(xG);
            dNdxTrial = trial.evaluateCartesianDerivatives(xG);
            nNode = size(dNdxTrial,2);
            nGaus = size(dNdxTrial,3);
            nElem = size(dNdxTrial,4);

            dofToDim = repmat(1:3,[1,nNode]);
            dofToNode = repmat(1:nNode,[1,3]);

            K = zeros(3*nNode,3*nNode,nGaus,nElem);
            for iDof = 1:24 % test dof
                iNode = dofToNode(iDof);
                iDim  = dofToDim(iDof);
                GradDeltaV = zeros(3,3, nGaus, nElem);
                GradDeltaV(iDim,:,:,:) = dNdxTest(:,iNode,:,:);

                res = zeros(3,3,nGaus,nElem);
                for a = 1:3
                    for b = 1:3
                        C = squeeze(Ctan(:,:,a,b,:,:));
                        res(a,b,:,:) = bsxfun(@(A,B) sum(A.*B, [1 2]), GradDeltaV,C);
                    end
                end

                for jDof = 1:24 % trial dof
                    jNode = dofToNode(jDof);
                    jDim  = dofToDim(jDof);

                    GradDeltaU = zeros(3,3, nGaus, nElem);
                    GradDeltaU(jDim,:,:,:) = dNdxTrial(:,jNode,:,:);
                    K(iDof,jDof,:,:) = bsxfun(@(A,B) sum(A.*B, [1 2]), res, GradDeltaU);

%                     for a = 1:3
%                         for b = 1:3
%                             C = squeeze(Ctan(a,b,:,:,:,:));
%                             res(a,b,:,:) = bsxfun(@(A,B) sum(A.*B, [1 2]), C,GradDeltaU);
%                         end
%                     end

%                     K(iDof,jDof,:,:) = bsxfun(@(A,B) sum(A.*B, [1 2]), GradDeltaV, res);
                end
            end
            K = K.*dV;
            K = squeeze(sum(K,3));

            s.fun    = []; % !!!
            assembler = AssemblerFun(s);
            hess = assembler.assemble(K, test, trial);
            % - Not symmetric
            % - Negative values at diagonal
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

            invFt = MatrixVectorizedInverter.computeInverse(Ft);
            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);

            piola = obj.mu*(F-invFt) + obj.lambda*log(jac).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun,xG)
            [F,I33] = obj.computeDeformationGradient(uFun, xG);
            Ft = permute(F, [2 1 3 4]);

            invF = MatrixVectorizedInverter.computeInverse(F);
            invFt = MatrixVectorizedInverter.computeInverse(Ft);

            jac(1,1,:,:)  = MatrixVectorizedInverter.computeDeterminant(F);
            logJac(1,1,:,:,:,:) = log(jac);

            Aneo = obj.lambda*obj.outerProduct(invFt, invFt) + ...
                obj.mu*obj.kron_topF(I33,I33) + ...
                (obj.mu-obj.lambda*logJac).*obj.kron_botF(invFt, invF);
        end

        function [F,I33] = computeDeformationGradient(obj, uFun, xG)
            nPoints  = size(xG,2);
            nElem = obj.mesh.nelem;
            nDimG = obj.mesh.ndim;
            nDimf = uFun.ndimf;

            GradUT = reshape(Grad(uFun).evaluate(xG),[nDimG,nDimf,nPoints, nElem]);
            GradU = permute(GradUT, [2 1 3 4]);
%             GradU = [0.0 0.0 0.0; -3.415063509461096 -0.24999999999999956 -0.4330127018922192; 0.9150635094610968 0.43301270189221924 -0.24999999999999994];

            I33 = zeros(size(GradU));
            I33(1,1,:,:) = 1;
            I33(2,2,:,:) = 1;
            I33(3,3,:,:) = 1;

            F = I33 + GradU;
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

        function C= kron_topF(obj, A,B)
            C = zeros([size(A,1), size(A,2), size(B)]); % to support 4th order tensors
            for i = 1:size(A,1)
                for k = 1:size(B,1)
                    C(i,:,k,:,:,:) = pagemtimes( A(i,k,:,:), B);
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