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
            % dtdt.evaluate(xV) 1 x ngaussxnelem
            % Expand(dtdt.evaluate(xV)) 1 x 1 x ngaussxnelem
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) %tangent

            % jump.setFValues([0,0.0505;0,0.0505]); %Tests

            d    = obj.computeDamage(jump);           % correct
            ddot = obj.computeDamageDerivative(jump); % correct
            unoZero = ConstantFunction.create([1,0],jump.mesh);
            zeroUno = ConstantFunction.create([0,1],jump.mesh);
            dtdt = obj.K * ((1-DP(d,unoZero)) - Expand(DP(ddot,unoZero).*DP(jump,unoZero),2));
            dndn = obj.K * ((1-DP(d,zeroUno)) - Expand(DP(ddot,zeroUno).*DP(jump,zeroUno),2));
            
            dtdt = dtdt.evaluate(xV);  %correct
            dndn = dndn.evaluate(xV); %correct
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            gradT =  [dtdt                    ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem), dndn                  ]; %size??
            % dtdt.evaluate(xV) 1 x ngaussxnelem
            % Expand(dtdt.evaluate(xV)) 1 x 1 x ngaussxnelem
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