classdef Grad < MathEntity
    
    properties (Access = private)
        F
    end
    
    methods (Access = public)
        
        function A = Grad(F)
            obj.init(F)
            A.ops = @(x) obj.evaluate(x);
            
        end
        
        function GradF = evaluate(obj,xV)
            dNdx = obj.F.evaluateCartesianDerivatives(xV);
            nDimf = obj.F.ndimf;
            nDimG = size(dNdx, 1);
            nNodeE = size(dNdx, 2);
            nPoints = size(dNdx, 3);
            nElem = size(dNdx, 4);
           
            fV = reshape(obj.F.fValues',[numel(obj.F.fValues) 1]);
            grad = zeros(nDimG, nDimf, nPoints, nElem);
            for iDimG = 1:nDimG
                for jDimf = 1:nDimf
                    for kNodeE = 1:nNodeE
                        dNdxIK = squeezeParticular(dNdx(iDimG, kNodeE,:,:),[1 2]);
                        iDofE = nDimf*(kNodeE-1)+jDimf;
                        dofs = obj.F.connec(:,iDofE);
                        fKJ = repmat(fV(dofs),[nPoints 1]);
                        gradIJ= dNdxIK.*fKJ;
                        grad(iDimG,jDimf,:,:) = squeezeParticular(grad(iDimG,jDimf,:,:),[1 2]) + gradIJ;
                    end
                end
            end
            GradF = grad;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,F)
            obj.F = F;
        end
        
    end
    
end