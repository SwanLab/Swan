classdef TractionBiliniarCoupled < handle

    properties (Access = private)
        jumpCrit
        fractureStrength %tau_o
        fractureToughness % g_c
    end

    properties (Access = private)
        K
        jumpFinal
    end
    
    methods (Access = public)
        
        function obj = TractionBiliniarCoupled(cParams)
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
            
            %arreglar parametres
            obj.fractureStrength  = cParams.fractureStrength;
            obj.fractureToughness = cParams.fractureToughness;
            obj.jumpCrit          = cParams.jumpCrit;
            obj.jumpFinal = 2*obj.fractureToughness/obj.fractureStrength;
            obj.K         = obj.fractureStrength/obj.jumpCrit;
        end

        function gradT = computeSecantGradientMatrix(obj, jump, xV) %secant

        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) %tangent

        end

        function d = computeDamage(obj,jump)
            jumpNorm = computeJumpNorm(jump);
            L = min(1, (obj.jumpFinal*(jumpNorm-obj.jumpCrit))/ ...
                (jumpNorm*(obj.jumpFinal-obj.jumpCrit)));
            d = max(0,L);
        end

        function jN = computeJumpNorm(obj,jump)
            unoZero = ConstantFunction.create([1,0],jump.mesh);
            zeroUno = ConstantFunction.create([0,1],jump.mesh);
            jN = (unoZero * jump) ^2 + (zeroUno * jump)^2;
            jN = sqrt(jN);
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