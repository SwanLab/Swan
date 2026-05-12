classdef TractionBiliniarCoupled < handle

    properties (Access = private)
        tau0 % fracture strength
        Gc   % fracture toughness
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
            jumpNorm = obj.computeJumpNorm(jump);
            d = obj.computeDamage(jumpNorm);
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
            obj.jumpCrit  = cParams.jumpCrit;
            obj.jumpFinal = cParams.jumpFinal;
            obj.K         = 1e8;
            % obj.jumpFinal = 2*obj.Gc/obj.tau0;
            % obj.jumpCrit  = obj.tau0/obj.K;
        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) % 2 x 2 x ngauss x nelem
            [dtdt,dtdn,dndt,dndn] = obj.computeTractionDerivatives(jump);
            dtdt=dtdt.evaluate(xV); dtdn=dtdn.evaluate(xV);
            dndt=dndt.evaluate(xV); dndn=dndn.evaluate(xV);
            ngauss = size(dtdt,2); nelem  = size(dtdt,3);
            gradT = [reshape(dtdt,1,1,ngauss,nelem),reshape(dtdn,1,1,ngauss,nelem);
                    reshape(dndt,1,1,ngauss,nelem),reshape(dndn,1,1,ngauss,nelem)];
        end

        function [dtdt,dtdn,dndt,dndn] = computeTractionDerivatives(obj, jump)
            isDamaging = obj.isJumpDamaging(jump);  % nDimJumpNorm (1) x ngauss x nelem
            jumpNorm = obj.computeJumpNorm(jump); % 1 x ngauss x nelem
            d    = obj.computeDamage(jumpNorm);             % 1 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump,jumpNorm,isDamaging);   % 2 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            ddot_t = DP(ddot,unoZero); ddot_n = DP(ddot,zeroUno);
            jumpT = DP(jump,unoZero); jumpN = DP(jump,zeroUno);
            dtdt = obj.K * ((1-d) - jumpT.*ddot_t); % 1 x ngauss x nelem
            dtdn = obj.K * (-jumpT.*ddot_n);
            dndt = obj.K * (-jumpN.*ddot_t);
            dndn = obj.K * ((1-d) - jumpN.*ddot_n);
        end

        function d = computeDamage(obj,jumpNorm) % 1 x ngauss x nelem   
            d = min((obj.jumpFinal*(jumpNorm-obj.jumpCrit))./ ...
                (jumpNorm*(obj.jumpFinal-obj.jumpCrit)),1);
            d = max(d,0);
            fprintf('d range: [%e , %e]\n',min(d.evaluate([-1,1]),[],'all'), max(d.evaluate([-1,1]),[],'all'));
        end

        function ddot = computeDamageDerivative(obj,jump,jumpNorm,isDamaging)
            alpha      = obj.jumpCrit * obj.jumpFinal / (obj.jumpCrit-obj.jumpFinal) ./jumpNorm^3;
            ddot       = alpha .* jump .* isDamaging;
        end

        function isDamaging = isJumpDamaging(obj,jump) 
            jumpNorm = obj.computeJumpNorm(jump)-1e-15;
            temp1 = jumpNorm - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - jumpNorm; % b - f
            isDamaging = temp1.*temp2 > 0;  % 1 x ngauss x nelem
        end       
    
        function jN = computeJumpNorm(obj,jump) % 1 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            jN = sqrt(DP(jump,unoZero,1,1).^2 + DP(jump,zeroUno,1,1).^2) + 1e-15;
        end
    end
end