classdef TractionBiliniarUncoupled < handle

    properties (Access = private)
        K
        j0
        jF
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
            d = (abs(jump)-obj.j0)./(obj.jF-obj.j0);
            d = max(min(d,1),0);
            % fprintf('d range: [%e , %e]\n', min(d.evaluate([-1,1]),[],'all'), max(d.evaluate([-1,1]),[],'all'));
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.K = cParams.Kcoh;
            obj.j0 = cParams.jumpCrit;
            obj.jF = cParams.jumpFinal;
        end

        function ddot =  computeDamageDerivative(obj,jump)
            % isDamaging = isJumpDamaging(obj,jump);
            d    =  obj.computeDamage(jump); % 1 x ngauss x nelem
            isNotFullyDamage = ( d < 1 );

            isDamaging = ((abs(jump)-obj.j0).*(abs(jump)-obj.jF)) < 0;

            ddot = (1./(obj.jF - obj.j0)).*isDamaging;
        end

    end
end