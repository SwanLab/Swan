classdef RHSintegrator_ShapeSymmDerivative < RHSintegrator

    properties (Access = private)
        test
        BdV
    end
    
    methods (Access = public)
        
        function obj = RHSintegrator_ShapeSymmDerivative(cParams)
            obj.init(cParams)
            obj.setQuadratureOrder(cParams);
            obj.createQuadrature();
            xV = obj.quadrature.posgp;
            dNdx = obj.test.evaluateCartesianDerivatives(xV);
            B    = obj.computeB(dNdx);  
            dV = obj.mesh.computeDvolume(obj.quadrature); 
            obj.BdV  = B.* reshape(dV, [1, 1, obj.mesh.nelem, size(xV,2)]);
        end
        
        function rhs = compute(obj, fun)
            rhsElem = obj.computeElementalRHS(fun);
            rhs = obj.assembleIntegrand(rhsElem);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.test = cParams.test;
        end
        %B = obj.B;
        function rhsC = computeElementalRHS(obj, fun)          
           %  [nDim, nNode, nGaus, nElem] = size(obj.dNdx);
           %  B = obj.B;%rand(nVoigt, nDofE, nElem, nGaus);
           %  fG = fun.evaluate(obj.quadrature.posgp);
           %  dV = obj.mesh.computeDvolume(obj.quadrature);            
           % 
           %  fG_reshaped = reshape(fG, size(fG, 1), nGaus, nElem);
           % fGI = permute(fG_reshaped, [1, 3, 2]);
           % fdv = fGI .* permute(dV, [3, 2, 1]);
           % fdv = reshape(fdv, [1, size(fdv, 1), nElem, nGaus]);
           % B_reshaped = permute(B, [1, 2, 3, 4]);
           % intI = pagemtimes(fdv, B_reshaped);
           % rhsC = sum(squeeze(intI), 3);

           fG = fun.evaluate(obj.quadrature.posgp);           
           fGP = permute(fG, [4, 1, 3, 2]);
           int = pagemtimes(fGP,obj.BdV);          
           rhsC = sum(squeeze(int),3);

           % fGP  = permute(fG, [4, 1, 3, 2]);
           % int  = pagemtimes(fGP, BdV);
           % rhsC = squeeze(int);

            % A = rand(3,4,5,3);
            % B = rand(3,4,1,3);
            % C = tensorprod(A,B,[1 2],[1 2])
          % norm(rhsC(:) - rhsC2(:))
        end

        function f = assembleIntegrand(obj, rhsElem)
            integrand = pagetranspose(rhsElem);
            connec = obj.test.getDofConnec();
            nDofs = max(max(connec));
            nDofElem  = size(connec,2);
            f = zeros(nDofs,1);
            for idof = 1:nDofElem
                int = integrand(:,idof);
                con = connec(:,idof);
                f = f + accumarray(con,int,[nDofs,1],@sum,0);
            end
        end

        function B = computeB(obj, dNdx)
            nGaus = obj.quadrature.ngaus;
            nNode = size(dNdx, 2);
            nElem = obj.mesh.nelem;
            nDim  = obj.mesh.ndim;
            nDimF = obj.test.ndimf;
            nDofE = nNode*nDimF;
            nVoigt = nDim * (nDim + 1) / 2;
            j = nDimF * reshape((1:nNode) - 1, 1, nNode, 1);
            B = zeros(nVoigt, nDofE, nElem, nGaus);
            d = permute(dNdx, [1, 2, 4, 3]);
            if nDim == 2
                B(1, j + 1, :, :) = d(1, :, :, :);
                B(2, j + 2, :, :) = d(2, :, :, :);
                B(3, j + 1, :, :) = d(2, :, :, :);
                B(3, j + 2, :, :) = d(1, :, :, :);
            elseif nDim == 3
                B(1, j + 1, :, :) = d(1, :, :, :);
                B(2, j + 2, :, :) = d(2, :, :, :);
                B(3, j + 3, :, :) = d(3, :, :, :);
                B(4, j + 1, :, :) = d(2, :, :, :);
                B(4, j + 2, :, :) = d(1, :, :, :);
                B(5, j + 1, :, :) = d(3, :, :, :);
                B(5, j + 3, :, :) = d(1, :, :, :);
                B(6, j + 2, :, :) = d(3, :, :, :);
                B(6, j + 3, :, :) = d(2, :, :, :);
            end
        end


    end
    
end