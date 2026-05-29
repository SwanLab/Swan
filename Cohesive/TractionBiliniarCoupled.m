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

        lambdaOld
        lambdaTrial
    end

    properties (Access = private)
        K
        jumpFinal
        jumpCrit




        trackLambda
        trackD
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

        function computeJumpNorm(obj,jump) % 1 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            jumpT = DP(jump,unoZero,1,1);
            jumpN = DP(jump,zeroUno,1,1);
            jumpNPos = 0.5*(jumpN + abs(jumpN));
            jN = sqrt(jumpT.^2 + jumpNPos.^2) + 1e-15;
            if isempty(obj.lambdaOld)
                obj.lambdaOld = zeros(size(jN.evaluate([-1,1])));
            end
            obj.lambdaTrial = max(jN,obj.lambdaOld);
        end

        function d = computeDamage(obj,jump)
            obj.computeJumpNorm(jump);
            obj.computeJumpLimits(jump);
            d = min((obj.jumpFinal .* (obj.lambdaTrial - obj.jumpCrit)./ (obj.lambdaTrial .* (obj.jumpFinal - obj.jumpCrit))),1);
            d = max(d,0);

            obj.trackD      = [obj.trackD;d.evaluate([-1,1])];
            obj.trackLambda = [obj.trackLambda;obj.lambdaTrial.evaluate([-1,1])];

            % fprintf('d range: [%e , %e]\n', min(d.evaluate([-1,1]),[],'all'), max(d.evaluate([-1,1]),[],'all'));
        end
 
        function updateLambdaOld(obj,lOld)
            obj.lambdaOld = obj.lambdaTrial;
        end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.K                = 3.2e13;
            obj.tau0Normal       = cParams.tau0Normal;
            obj.tau0Shear        = cParams.tau0Shear;
            obj.firstCritEnergy  = cParams.firstCritEnergy;
            obj.secondCritEnergy = cParams.secondCritEnergy;
            obj.eta              = cParams.eta;

            obj.jump0Normal = obj.tau0Normal / obj.K;
            obj.jump0Shear  = obj.tau0Shear / obj.K;

            obj.jumpFinalNormal = 2*obj.firstCritEnergy  / obj.tau0Normal;
            obj.jumpFinalShear  = 2*obj.secondCritEnergy / obj.tau0Shear;
            fprintf('jump0Normal = %e\n',obj.jump0Normal)
            fprintf('jumpFinalNormal = %e\n',obj.jumpFinalNormal)


            obj.trackD = zeros(1,2);
            obj.trackLambda = obj.trackD;

        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) % 2 x 2 x ngauss x nelem
            obj.computeJumpNorm(jump);
            [dtdt,dtdn,dndt,dndn] = obj.computeTractionDerivatives(jump);
            dtdt=dtdt.evaluate(xV); dtdn=dtdn.evaluate(xV);
            dndt=dndt.evaluate(xV); dndn=dndn.evaluate(xV);
            ngauss = size(dtdt,2); nelem  = size(dtdt,3);
            gradT = [reshape(dtdt,1,1,ngauss,nelem),reshape(dtdn,1,1,ngauss,nelem);
                    reshape(dndt,1,1,ngauss,nelem),reshape(dndn,1,1,ngauss,nelem)];
        end

        function [dtdt,dtdn,dndt,dndn] = computeTractionDerivatives(obj, jump)
            obj.computeJumpLimits(jump);
            d    = obj.computeDamage(jump);             % 1 x ngauss x nelem
            ddot = obj.computeDamageDerivative(jump);   % 2 x ngauss x nelem
            unoZero = ConstantFunction.create([1;0],jump.mesh);
            zeroUno = ConstantFunction.create([0;1],jump.mesh);
            ddot_t = DP(ddot,unoZero); ddot_n = DP(ddot,zeroUno);
            jumpT = DP(jump,unoZero); jumpN = DP(jump,zeroUno);
            dtdt = obj.K * ((1-d) - jumpT.*ddot_t); % 1 x ngauss x nelem
            dtdn = obj.K * (-jumpT.*ddot_n);
            dndt = obj.K * (-jumpN.*ddot_t);
            dndn = obj.K * ((1-d) - jumpN.*ddot_n);
        end

        function ddot = computeDamageDerivative(obj,jump)
            isDerivativeZero = obj.checkIsDerivativeZero(jump);  % nDimJumpNorm (1) x ngauss x nelem
            obj.computeJumpNorm(jump);
            alpha      = obj.jumpCrit .* obj.jumpFinal ./ (obj.jumpFinal-obj.jumpCrit) ./obj.lambdaTrial.^3;
            ddot       = alpha .* jump .* isDerivativeZero;
        end

        function isDerivativeZero = checkIsDerivativeZero(obj,jump) 
            obj.computeJumpNorm(jump);
            temp1 = obj.lambdaTrial - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - obj.lambdaTrial; % b - f
            isDerivativeZero = temp1.*temp2 > 0;  % 1 x ngauss x nelem
        end       

        function B = computeMixedModeRatio(obj,jump,isJump) % jumpNorm = 1xngaussxnelem | jumpshear = 1xngaussxnelem | B = 1xngaussxnelem 
            unoZero   = ConstantFunction.create([1;0],jump.mesh);
            jumpShear = DP(jump.',unoZero); 
            B         = isJump .* (jumpShear./obj.lambdaTrial).^2;
        end

        function computeJumpLimits(obj,jump)
            isJump = obj.lambdaTrial > 1e-12;
            B = obj.computeMixedModeRatio(jump,isJump);
                B = 0;      % Impose mode 1 check
            obj.jumpCrit  = sqrt(obj.jump0Normal.^2 + (obj.jump0Shear.^2 - obj.jump0Normal.^2).*B.^obj.eta);
            obj.jumpFinal = (obj.jump0Normal*obj.jumpFinalNormal + (obj.jump0Shear*obj.jumpFinalShear - obj.jump0Normal*obj.jumpFinalNormal).*B.^obj.eta)./obj.jumpCrit;
            % 
            % obj.jumpCrit  = 1.25e-7; 
            % obj.jumpFinal = 0.025e-3;
        end
    end
end