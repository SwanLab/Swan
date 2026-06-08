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
            
            % Printing ========================
                lambda = obj.computeJumpNorm(jump);
                [j0, jF] = obj.computeJumpLimits(jump,lambda);
                ddot = j0 .* jF ./ ((jF-j0) .*lambda.^3);
                d    =  obj.computeDamage(jump); % 1 x ngauss x nelem
                fprintf('lambda = %.3e\n', lambda.evaluate([-1,1]));
                fprintf('j0     = %.3e\n', j0.evaluate([-1,1]));
                fprintf('jF     = %.3e\n', jF.evaluate([-1,1]));
                fprintf('d      = %.3e\n', max(d.evaluate([-1,1])));
            % =================================


            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            dtSec = obj.computeDerivativeSecant(jump);
            dtTan = dtSec -  obj.K * ddot.*kronProd(jump,jump,[1 2]);

        end 

        function d = computeDamage(obj,jump)
            lambda = obj.computeJumpNorm(jump);
            [j0, jF] = obj.computeJumpLimits(jump,lambda);
            d = (jF .* (lambda - j0)) ./ (lambda .* (jF - j0)) ;
            d = max(min(d,1),obj.dOld);
            obj.dOld = d.evaluate([-1,1]);
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

        % function dd = computeDamageDerivative(obj,jump)
        % 
        %     lambda = obj.computeJumpNorm(jump);
        %     [j0,jF] = obj.computeJumpLimits(jump,lambda);
        % 
        %     unoZero = ConstantFunction.create([1;0],jump.mesh);
        %     zeroUno = ConstantFunction.create([0;1],jump.mesh);
        % 
        %     Dt = DP(jump.',unoZero);   % shear
        %     Dn = DP(jump.',zeroUno);   % normal
        % 
        %     Dnp = max(Dn,0);
        % 
        %     dl_dDt = Dt./lambda;
        %     dl_dDn = Dnp./lambda;
        % 
        %     B = (Dt./lambda).^2;
        % 
        %     dB_dDt = 2*Dt./lambda.^2 ...
        %            - 2*Dt.^3./lambda.^4;
        % 
        %     dB_dDn = -2*Dt.^2 .* Dnp ./ lambda.^4;
        % 
        %     A0 = obj.jump0Shear.^2 - obj.jump0Normal.^2;
        % 
        %     dj0_dB = ...
        %         (A0 .* obj.eta .* B.^(obj.eta-1)) ...
        %         ./ (2*j0);
        % 
        %     dj0_dDt = dj0_dB .* dB_dDt;
        %     dj0_dDn = dj0_dB .* dB_dDn;
        % 
        %     C0 = obj.jump0Normal*obj.jumpFinalNormal;
        % 
        %     C1 = obj.jump0Shear*obj.jumpFinalShear ...
        %        - obj.jump0Normal*obj.jumpFinalNormal;
        % 
        %     N = C0 + C1.*B.^obj.eta;
        % 
        %     dN_dB = C1 .* obj.eta .* B.^(obj.eta-1);
        % 
        %     djF_dB = ...
        %         (dN_dB .* j0 - N .* dj0_dB) ...
        %         ./ (j0.^2);
        % 
        %     djF_dDt = djF_dB .* dB_dDt;
        %     djF_dDn = djF_dB .* dB_dDn;
        % 
        %     dd_dlambda = ...
        %         (j0 .* jF) ...
        %         ./ ((jF-j0).*lambda.^2);
        % 
        %     dd_dj0 = ...
        %         jF .* (lambda-jF) ...
        %         ./ (lambda .* (jF-j0).^2);
        % 
        %     dd_djF = ...
        %         j0 .* (j0-lambda) ...
        %         ./ (lambda .* (jF-j0).^2);
        % 
        %     dd_dDt = ...
        %         dd_dlambda .* dl_dDt + ...
        %         dd_dj0     .* dj0_dDt + ...
        %         dd_djF     .* djF_dDt;
        % 
        %     dd_dDn = ...
        %         dd_dlambda .* dl_dDn + ...
        %         dd_dj0     .* dj0_dDn + ...
        %         dd_djF     .* djF_dDn;
        % 
        %     dd = [dd_dDt,dd_dDn];
        %     d = obj.computeDamage(jump);
        %     dd = dd .* (d < 1);
        % end

        function lambda = computeJumpNorm(~,jump) % 1 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            jumpT = DP(jump,unoZero,1,1);
            jumpN = DP(jump,zeroUno,1,1);
            lambda = sqrt(jumpT.^2 + macaulay(jumpN).^2) + 1e-10;           
        end

        function [jump0, jumpFinal] = computeJumpLimits(obj,jump,lambda)
            B = obj.computeMixedModeRatio(jump,lambda);
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