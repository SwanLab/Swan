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
            obj.Gc    = cParams.fractureToughness;
            obj.jumpCrit          = cParams.jumpCrit;
            obj.K         = 1e8;
            obj.jumpFinal = 2*obj.Gc/obj.tau0;
            obj.jumpCrit  = obj.tau0/obj.K;
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) % 2 x 2 x ngauss x nelem
            d    = obj.computeDamage(jump);             % 2 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            unoZero = ConstantFunction.create([1,0],jump.mesh);
            zeroUno = ConstantFunction.create([0,1],jump.mesh);
            ddot_t = DP(ddot,unoZero);
            ddot_n = DP(ddot,zeroUno);
            jumpT = DP(jump,unoZero);
            jumpN = DP(jump,zeroUno);
 
            dtdt = obj.K * ((1-d) - jumpT.*ddot_t); % mida malament (hauria de ser 1 x 1 x 2
            dtdn = obj.K * (-jumpT.*ddot_n);
            dndt = obj.K * (-jumpN.*ddot_t);
            dndn = obj.K * ((1-d) - jumpN.*ddot_n);

            dtdt=dtdt.evaluate(xV); dtdn=dtdn.evaluate(xV);
            dndt=dndt.evaluate(xV); dndn=dndn.evaluate(xV);

            ngauss = size(dtdt,3); 
            nelem  = size(dtdt,4);

            tmp = zeros(2,ngauss,nelem);
            dtdt=tmp;dtdn=tmp;dndt = tmp;dndn=tmp;
            
            % 1 x 1 x ngauss x nelem
            dtdt = reshape(dtdt,1,1,ngauss,nelem);
            dtdn = reshape(dtdn,1,1,ngauss,nelem);
            dndt = reshape(dndt,1,1,ngauss,nelem);
            dndn = reshape(dndn,1,1,ngauss,nelem);
            
            gradT = [dtdt,dtdn; dndt,dndn];

        end

        function d = computeDamage(obj,jump) % 1 x ngauss x nelem
            jumpNorm = obj.computeJumpNorm(jump);
            L = min((obj.jumpFinal*(jumpNorm-obj.jumpCrit))./ ...
                (jumpNorm*(obj.jumpFinal-obj.jumpCrit)),1);
            d = max(L,0); % 1 x ngauss x nelem
        end

        function jN = computeJumpNorm(obj,jump) % 1 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            jN = sqrt(DP(jump,unoZero,1,1).^2 + DP(jump,zeroUno,1,1).^2);
        end

        function ddot =  computeDamageDerivative(obj,jump)
            isDamaging = obj.isJumpDamaging(jump);
            jumpNorm = obj.computeJumpNorm(jump);
            alpha    = -obj.jumpCrit * obj.jumpFinal / (obj.jumpCrit-obj.jumpFinal) ./ max(jumpNorm^3,1e-8);
            ddot     = alpha .* jump;
            ddot     = ddot.*isDamaging;
        end

        function isDamaging = isJumpDamaging(obj,jump) % comprovar!!
            temp1 = jump - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - jump; % b - f
            isDamaging = temp1.*temp2 > 0; 
        end        
    end
end