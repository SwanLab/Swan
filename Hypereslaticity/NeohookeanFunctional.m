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
        end

        function val = compute(obj, uFun)      
            nDim = obj.mesh.ndim;
            [~,F] = obj.computeDeformationGradient(uFun);
            C = F'*F;
            trC = trace(C);
            jac = Det(F);
            s.quadType = 2;
            s.mesh = obj.mesh;           
            int = IntegratorFunction(s);
            f = obj.mu./2.*(trC - nDim) - obj.mu.*log(jac) + obj.lambda./2.*(log(jac)).^2;
            val = int.compute(f);
        end
        
        function Fint = computeGradient(obj, uFun)
%             piola = PK1.evaluate(xG); 
%             dofToDim = repmat(1:nDimf,[1,nNode]);
%             dofToNode = repmat(1:nNode, nDimf, 1);
%             dofToNode = dofToNode(:);
% 
%             fint = zeros(nDof,1,nGaus,nElem);
%             % fint2 = zeros(nDof,1,nGaus,nElem);
% 
% 
%             % for iNode = 1:nNode
%             %     for iDim = 1:nDim
%             %         for iField = 1:nDimf
%             %             dNdx_ij = dNdxTest(iDim, iNode, :, :);
%             %             iDof = iField + nDimf*(iNode-1);
%             %             Pik = piola(iField,iDim,:,:);
%             %             fint2(iDof,1,:,:) = fint2(iDof,1,:,:) + Pik.*dNdx_ij;
%             %         end
%             %     end
%             % end
% 
%             % GradDeltaV is not always compatible (see BCs), but we dont
%             % worry about it since we reduce the matrix later on    
%             for iDof = 1:nDof
%                 iNode = dofToNode(iDof);
%                 iDim  = dofToDim(iDof);
%                 deltav = zeros(nNode,nDimf, nGaus, nElem);
%                 deltav(iNode,iDim,:,:) = 1;
%                 GradDeltaV = pagemtimes(dNdxTest,deltav);
%                 fint(iDof, :,:,:) = squeeze(bsxfun(@(A,B) sum(A.*B, [1 2]), pagetranspose(piola),GradDeltaV));
%             end
%             % err = norm(fint(:)-fint2(:))/norm(fint(:))
%             % fint = fint2;
%             fint = fint.*dV;
%             fint = squeeze(sum(fint,3));
%             Fint = obj.assembleIntegrand(fint,test);
            
            PK1 = obj.computeFirstPiola(uFun);
            s.mesh = obj.mesh;
            s.quadratureOrder = 3;
            s.type = 'ShapeDerivativeTensor';
            rhs = RHSIntegrator.create(s);
            Fint = rhs.compute(PK1,uFun); 
        end

        function f = assembleIntegrand(obj, rhsElem, test)
            integrand = pagetranspose(rhsElem);
            connec = test.getDofConnec();
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
%             % This is the LINEALIZED hessian (Holzapfel, 401)
%             % See  Holzapfel, 396
            Aneofun = obj.computeTangentConstitutive(uFun);
            s.quadratureOrder = 3;
            s.test  = uFun;
            s.trial = uFun;
            s.mesh  = obj.mesh;
            s.type ='StiffnessFiniteStrain';
            lhs = LHSIntegrator.create(s);
            hess = lhs.compute(Aneofun);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.lambda = cParams.material.lambda;
            obj.mu     = cParams.material.mu;
            obj.mesh   = cParams.mesh;
        end

        function PK1 = computeFirstPiola(obj,uFun,xG)
            [~,F] = obj.computeDeformationGradient(uFun);
             invFt = Inv(F');
             jac = Det(F);
             PK1 = obj.mu.*(F-invFt) + obj.lambda.*log(jac).*invFt;
        end

        function Aneo = computeTangentConstitutive(obj,uFun)
            [I33,F] = obj.computeDeformationGradient(uFun);
            invF = Inv(F);
            jac = Det(F);
            Aneo = obj.lambda.*OP(invF', invF') + ...
                obj.mu.*kronTop(I33,I33) + ...
                Expand((obj.mu-obj.lambda.*log(jac)),4).*kronBot(invF', invF);
        end

        function [I33,F] = computeDeformationGradient(obj, uFun)
            I33 = Identity(uFun);
            F = I33 + Grad(uFun)';
        end     
    end
    
end