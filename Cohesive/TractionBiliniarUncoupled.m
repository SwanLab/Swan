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

        function gradT = computeTangentGradientMatrix(obj, jump, xV) %tangent
            d = obj.computeDamage(jump);
            dddj = obj.computeDamageDerivative(jump);
            dTdT = obj.K * ((1-DP(d,[1,0])) - DP(dddj,[1,0]) * DP(jump,[1 0]));%igual s'ha d'arreglar això (potser repmat)
            dNdN = obj.K * ((1-DP(d,[0,1])) - DP(dddj,[0,1]) * DP(jump,[1 0]));
            dtdt = Expand(dTdT,2).evaluate(xV); 
            dndn = Expand(dNdN,2).evaluate(xV);
            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);
            gradT =  [dtdt                    ,zeros(1,1,ngauss,nelem);
                      zeros(1,1,ngauss,nelem), dndn                  ];
            % dtdt.evaluate(xV) 1 x ngaussxnelem
            % Expand(dtdt.evaluate(xV)) 1 x 1 x ngaussxnelem
        end

        function d = computeDamage(obj,jump)
            isOverLimit = isJumpOverLimit(obj,jump);
            isDamaging = isJumpDamaging(obj,jump);
            d =0+1*isOverLimit+(jump-obj.jumpCrit)/(obj.jumpFinal-obj.jumpCrit)*isDamaging;
        end


        function dd =  computeDamageDerivative(obj,jump)
            isDamaging = isJumpDamaging(obj,jump);
            dd = (1/(obj.jumpFinal - obj.jumpCrit)) * isDamaging;
        end


        function isOverLimit = isJumpOverLimit(obj,jump)
            isOverLimit = (jump > obj.jumpFinal);
        end

        function isUnderLimit = isJumpUnderLimit(obj,jump)
            isUnderLimit = (jump < obj.jumpCrit);
        end

        function isDamaging = isJumpDamaging(obj,jump) % comprovar!!
            tempJumpCrit = ConstantFunction.create(obj.jumpCrit,jump.mesh); 
            temp1 = jump - tempJumpCrit; % f - a
            tempJumpFinal = ConstantFunction.create(obj.jumpFinal,jump.mesh);
            temp2 = tempJumpFinal - jump; % b - f
            isDamaging = temp1*temp2 > 0; % com evaluo isDamaging?
        end
        

        end
end