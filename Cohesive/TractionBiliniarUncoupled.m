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

            % unoZero = ConstantFunction.create([1,0],jump.mesh);
            % zeroUno = ConstantFunction.create([0,1],jump.mesh);
            % dtdt = obj.K * (1 - DP(d,unoZero)); 
            % dndn = obj.K * (1 - DP(d,zeroUno));
            % dtdt = dtdt.evaluate([-1,1]);
            % dndn = dndn.evaluate([-1,1]);
            % ngauss = size(dtdt,3); 
            % nelem  = size(dtdt,4);
            % gradT =  [dtdt                   ,zeros(1,1,ngauss,nelem);
            %           zeros(1,1,ngauss,nelem),dndn                  ];
        end

        function dtTan = computeDerivativeTangent(obj, jump) %2 x 2 x ngauss x nelem
            d    = obj.computeDamage(jump);             % 2 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            dtSec = (1-d).*obj.K.*eye(2);
            dtTan = dtSec - ddot.*squeezeParticular(kronProd(Expand(jump,2),Expand(jump,2),[1 2 3 4]),[1 2 3 4]);
        
            % unoZero = ConstantFunction.create([1,0],jump.mesh);
            % zeroUno = ConstantFunction.create([0,1],jump.mesh);
            % dtdt = obj.K * Expand((1-DP(d,unoZero)),2) - Expand(DP(ddot,unoZero).*DP(jump,unoZero),2); 
            % dndn = obj.K * Expand((1-DP(d,zeroUno)),2) - Expand(DP(ddot,zeroUno).*DP(jump,zeroUno),2);
            % dtdt = dtdt.evaluate([-1,1]);
            % dndn = dndn.evaluate([-1,1]);
            % ngauss = size(dtdt,3); 
            % nelem  = size(dtdt,4);
            % dtdt = reshape(dtdt,1,1,ngauss,nelem); % 1 x 1 x ngauss x nelem
            % dndn = reshape(dndn,1,1,ngauss,nelem); % 1 x 1 x ngauss x nelem
            % %2 x 2 x ngauss x nelem
            % gradT =  [dtdt                    ,zeros(1,1,ngauss,nelem);
            %           zeros(1,1,ngauss,nelem), dndn                  ]; 
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