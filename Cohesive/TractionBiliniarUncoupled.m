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
            s.operation = @(xV) obj.computeGradientMatrix(jump,xV);
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

        function d = computeDamage(obj,jump)
            isOverLimit = isJumpOverLimit(obj,jump);
            isUnderLimit = isJumpUnderLimit(obj,jump);
            isDamaging = isJumpDamaging(obj,jump);
            d = 1* isOverLimit+(jump-obj.jumpCrit)/(obj.jumpFinal - obj.jumpCrit)*isDamaging+isUnderLimit*0;
        end

        function gradT = computeGradientMatrix(obj, jump, xV) %secant
            d = obj.computeDamage(jump);
            dTdT = obj.K * (1-DP(d,[1,0]));%igual s'ha d'arreglar això (potser repmat)
            dNdN = obj.K * (1-DP(d,[0,1]));
            dtdt = Expand(dTdT,2).evaluate(xV); 
            dndn = Expand(dNdN,2).evaluate(xV);
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            gradT =  [dtdt                   ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem),dndn                  ];
            % dtdt.evaluate(xV) 1 x ngaussxnelem
            % Expand(dtdt.evaluate(xV)) 1 x 1 x ngaussxnelem
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) %secant
            d = obj.computeDamage(jump);
            dTdT = obj.K * (1-DP(d,[1,0]));%igual s'ha d'arreglar això (potser repmat)
            dNdN = obj.K * (1-DP(d,[0,1]));
            dtdt = Expand(dTdT,2).evaluate(xV); 
            dndn = Expand(dNdN,2).evaluate(xV);
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            gradT =  [dtdt                   ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem),dndn                  ];
            % dtdt.evaluate(xV) 1 x ngaussxnelem
            % Expand(dtdt.evaluate(xV)) 1 x 1 x ngaussxnelem
        end

        function isOverLimit = isJumpOverLimit(obj,jump)
            isOverLimit = (jump > obj.jumpFinal);
        end

        function isUnderLimit = isJumpUnderLimit(obj,jump)
            isUnderLimit = (jump < obj.jumpCrit);
        end

        function isDamaging = isJumpDamaging(obj,jump)            
            % tempJumpCrit =
            % ConstantFunction.create(obj.jumpCrit,jump.mesh); 
            % he de trobar una solucio per això
            % temp = jump - tempJumpCrit;
            % isDamaging = temp > 0;
        end
    end
end