classdef TractionBiliniarCoupled < handle

    properties (Access = private)
        tau0 % fracture strength
        Gc % fracture toughness

        initialFracture 
    end

    properties (Access = private)
        K
        jumpFinal
        jumpCrit
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
            obj.tau0  = cParams.fractureStrength;
            obj.Gc = cParams.fractureToughness;
            obj.jumpCrit          = cParams.jumpCrit;
            obj.K         = 1e8;
            obj.jumpFinal = 2*obj.Gc/obj.tau0;
            obj.jumpCrit  = obj.tau0/obj.K;
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) % 2 x 2 x ngauss x nelem tangent!
            
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

        end

        function isDamaging = isJumpDamaging(obj,jump) % comprovar!!
            temp1 = jump - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - jump; % b - f
            isDamaging = temp1.*temp2 > 0; 
        end

        
    end
end