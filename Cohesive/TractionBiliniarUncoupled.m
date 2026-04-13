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

        function dt = computeDerivative(obj,jump)
            s.operation = @(xV) obj.computeTangentGradientMatrix(jump,xV);
            s.ndimf = 4;
            s.mesh = jump.mesh;
            dt = DomainFunction(s);
        end
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.K = cParams.K;
            obj.jumpCrit  = cParams.jumpCrit;
            obj.jumpFinal = cParams.jumpFinal;
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) %2 x 2 x ngauss x nelem
            d    = obj.computeDamage(jump);             % 2 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            unoZero = ConstantFunction.create([1,0],jump.mesh);
            zeroUno = ConstantFunction.create([0,1],jump.mesh);
            dtdt = obj.K * Expand((1-DP(d,unoZero)),2) - Expand(DP(ddot,unoZero).*DP(jump,unoZero),2); 
            dndn = obj.K * Expand((1-DP(d,zeroUno)),2) - Expand(DP(ddot,zeroUno).*DP(jump,zeroUno),2);
            dtdt = dtdt.evaluate(xV);
            dndn = dndn.evaluate(xV);
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            dtdt = reshape(dtdt,1,1,ngauss,nelem); % 1 x 1 x ngauss x nelem
            dndn = reshape(dndn,1,1,ngauss,nelem); % 1 x 1 x ngauss x nelem
            %2 x 2 x ngauss x nelem
            gradT =  [dtdt                    ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem), dndn                  ]; 
        end

        function gradT = computeSecantGradientMatrix(obj, jump, xV) %secant
            d = obj.computeDamage(jump);
            unoZero = ConstantFunction.create([1,0],jump.mesh);
            zeroUno = ConstantFunction.create([0,1],jump.mesh);
            dTdT = obj.K * (1 - DP(d,unoZero)); 
            dNdN = obj.K * (1 - DP(d,zeroUno));
            dtdt = Expand(dTdT,2).evaluate(xV); 
            dndn = Expand(dNdN,2).evaluate(xV);
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            gradT =  [dtdt                   ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem),dndn                  ];
        end

        function d = computeDamage(obj,jump)
            d = (jump-obj.jumpCrit)./(obj.jumpFinal-obj.jumpCrit);
            d = max(min(d,1),0);
        end

        function ddot =  computeDamageDerivative(obj,jump)
            isDamaging = isJumpDamaging(obj,jump);
            ddot = (1./(obj.jumpFinal - obj.jumpCrit)).*isDamaging;
        end

        function isDamaging = isJumpDamaging(obj,jump) % comprovar!!
            temp1 = jump - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - jump; % b - f
            isDamaging = temp1.*temp2 > 0; 
        end
    end
end