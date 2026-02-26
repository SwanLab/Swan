classdef PreconditionerBlockDiagonal < handle
    
    properties (Access = private)
        LHSinv
    end
    
    properties (Access = private)        
        LHS
        dimension
        row
        sizeLHS
    end
    
    methods (Access = public)
        
        function obj = PreconditionerBlockDiagonal(cParams)
            obj.init(cParams)
            obj.computeBlockDiagonalInverseMatrix()
        end
        
        function z = apply(obj,r)            
            Pl = @(r) obj.LHSinv*r;
            z  = Pl(r);
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.LHS = cParams.LHS;
            obj.sizeLHS = size(obj.LHS,1);
            dim = cParams.dimension;
            obj.dimension = obj.closestDivisorMultiple(obj.sizeLHS,dim);
        end
        
        function computeBlockDiagonalInverseMatrix(obj)
            bLHS        = full(obj.LHS(1:obj.dimension,1:obj.dimension));
            bLHSinv     = pinv(bLHS);
            nrep        = obj.sizeLHS/obj.dimension;
            obj.LHSinv  = kron(speye(nrep), bLHSinv); 
%             obj.LHSinv  = repmat(bLHSinv,nrep);
        end

        function closest = closestDivisorMultiple(obj, num, target)
            candidates = 1:num;
            divisors = candidates(mod(num, candidates) == 0);

            % Find the divisor closest to the target
            [~, idx] = min(abs(divisors - target));
            closest = divisors(idx);
        end

        
    end
    
end