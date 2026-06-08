classdef TractionBiliniarUncoupled < handle

    properties (Access = private)
        K
        jumpCrit
        jumpFinal
    end
    
    methods (Access = public)
        
        function obj = TractionBiliniarUncoupled(cParams)
            obj.init(cParams)            
        end

        function t = computeFunction(obj,jump)
            d = obj.computeDamage(jump);
            t = obj.K * (1-d).*jump;
        end

        function dtSec = computeDerivativeSecant(obj, jump) %secant
            d = obj.computeDamage(jump);
            dtSec = (1-d).*obj.K.*eye(2);

        end

        function dtTan = computeDerivativeTangent(obj, jump) %2 x 2 x ngauss x nelem
            d    = obj.computeDamage(jump);             % 2 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            dtSec = (1-d).*obj.K.*eye(2);
            dtTan = dtSec - obj.K * ddot.*kronProd(jump,jump,[1 2]);

        end

        function d = computeDamage(obj,jump)
            d = (jump-obj.jumpCrit)./(obj.jumpFinal-obj.jumpCrit);
            d = max(min(d,1),0);
            % fprintf('d range: [%e , %e]\n', min(d.evaluate([-1,1]),[],'all'), max(d.evaluate([-1,1]),[],'all'));
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.K = cParams.Kcoh;
            obj.jumpCrit  = cParams.jumpCrit;
            obj.jumpFinal = cParams.jumpFinal;
        end

        function ddot =  computeDamageDerivative(obj,jump)
            % isDamaging = isJumpDamaging(obj,jump);
            d    =  obj.computeDamage(jump); % 1 x ngauss x nelem
            isNotFullyDamage = ( d < 1 );
            ddot = (1./(obj.jumpFinal - obj.jumpCrit)).*isNotFullyDamage;
        end

    end
end