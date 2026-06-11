classdef TractionBiliniarCoupled < handle

    properties (Access = private)
        tau0Normal
        tau0Shear
        firstCritEnergy
        secondCritEnergy
        eta
        jump0Normal
        jump0Shear
        jumpFinalNormal
        jumpFinalShear
    end

    properties (Access = private)
        K
        dOld
    end
    
    methods (Access = public)
        
        function obj = TractionBiliniarCoupled(cParams)
            obj.init(cParams)            
        end

        function t = computeFunction(obj,jump)
            d = obj.computeDamage(jump);
            t = obj.K * (1-d).*jump;
        end

        function dtSec = computeDerivativeSecant(obj,jump)
            d    =  obj.computeDamage(jump);              % 1 x ngauss x nelem
            dtSec = (1-d).*obj.K.*eye(2);
        end

        function dtTan = computeDerivativeTangent(obj,jump)
            
            % % Printing ========================
            %     lambda = obj.computeJumpNorm(jump);
            %     [j0, jF] = obj.computeJumpLimits(jump,lambda);
            %     ddot = j0 .* jF ./ ((jF-j0) .*lambda.^3);
            %     d    =  obj.computeDamage(jump); % 1 x ngauss x nelem
            %     fprintf('lambda = %.3e\n', lambda.evaluate([-1,1]));
            %     fprintf('j0     = %.3e\n', j0.evaluate([-1,1]));
            %     fprintf('jF     = %.3e\n', jF.evaluate([-1,1]));
            %     fprintf('d      = %.3e\n', max(d.evaluate([-1,1])));
            % % =================================

            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            dtSec = obj.computeDerivativeSecant(jump);
            dtTan = dtSec -  obj.K * ddot.*kronProd(jump,jump,[1 2]) + 1e-15;
            % dtTan = dtSec;
        end 

        function d = computeDamage(obj,jump)
            lambda = obj.computeJumpNorm(jump);
            [j0, jF] = obj.computeJumpLimits(jump,lambda);
            d = (jF .* (lambda - j0)) ./ (lambda .* (jF - j0)) ;
            d = max(min(d,1),obj.dOld);
        end

        function updateDamageOld(obj,jump)
            obj.dOld = obj.computeDamage(jump);
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.K                = cParams.Kcoh;
            obj.tau0Normal       = cParams.tau0Normal;
            obj.tau0Shear        = cParams.tau0Shear;
            obj.firstCritEnergy  = cParams.firstCritEnergy;
            obj.secondCritEnergy = cParams.secondCritEnergy;
            obj.eta              = cParams.eta;

            obj.jump0Normal = obj.tau0Normal / obj.K;
            obj.jump0Shear  = obj.tau0Shear / obj.K;

            obj.jumpFinalNormal = 2*obj.firstCritEnergy  / obj.tau0Normal;
            obj.jumpFinalShear  = 2*obj.secondCritEnergy / obj.tau0Shear;

            obj.dOld = 0;
        end

        function ddot = computeDamageDerivative(obj,jump)
            lambda = obj.computeJumpNorm(jump);
            [j0, jF] = obj.computeJumpLimits(jump,lambda);
            ddot = j0 .* jF ./ ((jF-j0) .*lambda.^3);

            isDamaging = ((lambda-j0).*(lambda-jF)) < 0;
            ddot = ddot .* isDamaging;
        end

        function lambda = computeJumpNorm(~,jump) % 1 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            jumpT = DP(jump,unoZero,1,1);
            jumpN = DP(jump,zeroUno,1,1);
            lambda = sqrt(jumpT.^2 + jumpN.^2) + 1e-10;           
        end

        function [jump0, jumpFinal] = computeJumpLimits(obj,jump,lambda)
            B = obj.computeMixedModeRatio(jump,lambda) ;
            jump0  = sqrt(obj.jump0Normal.^2 + (obj.jump0Shear.^2 - obj.jump0Normal.^2).*B.^obj.eta);
            jumpFinal = (obj.jump0Normal*obj.jumpFinalNormal + ...
                (obj.jump0Shear*obj.jumpFinalShear - obj.jump0Normal*obj.jumpFinalNormal).*B.^obj.eta)./jump0;
        end
            % jumpNorm = 1xngaussxnelem | jumpshear = 1xngaussxnelem | B = 1xngaussxnelem 
        function B = computeMixedModeRatio(~,jump,lambda) 
            unoZero   = ConstantFunction.create([1;0],jump.mesh);
            jumpShear = DP(jump.',unoZero);
            B         = (jumpShear./lambda).^2;
        end

    end
end