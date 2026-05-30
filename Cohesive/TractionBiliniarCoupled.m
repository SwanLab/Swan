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

            % obj.trackD      = [obj.trackD;d.evaluate([-1,1])];
            % obj.trackLambda = [obj.trackLambda;obj.lambdaTrial.evaluate([-1,1])];

            % fprintf('d range: [%e , %e]\n', min(d.evaluate([-1,1]),[],'all'), max(d.evaluate([-1,1]),[],'all'));
        end
 
        function updateLambdaOld(obj,lOld)
            obj.lambdaOld = obj.lambdaTrial;
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

            obj.trackD = zeros(1,2);
            obj.trackLambda = obj.trackD;

        end

        function gradT = computeTangentGradientMatrix(obj, jump, xV) % 2 x 2 x ngauss x nelem
            obj.computeJumpNorm(jump);
            obj.computeJumpLimits(jump);
            alpha = obj.computeDamageDerivative(jump,xV);   % 2 x ngauss x nelem
            d    = obj.computeDamage(jump);             % 1 x ngauss x nelem
            ngauss = size(d.evaluate([-1,1]),2);
            I = repmat(eye(2),1,1,ngauss,jump.mesh.nelem); % 2x2xngaussxnelem
            J = jump.evaluate(xV); % 2 x ngauss x nelem
            JJ = reshape(J,2,1,ngauss,jump.mesh.nelem).*reshape(J,1,2,ngauss,jump.mesh.nelem);
            gradT = (1-d) .* obj.K .* I +  alpha .* JJ;
            gradT = gradT.evaluate(xV);
            % gradT = (1-d)K*I - alpha * kronProd(jump,jump,[1 2 3 4])
            % (segurament amb squeeze al kron)
            % kronProd(sigBar,sigBar,[1 2 3 4]); hauria de sortir
            % 2x1x2xngaussxnelem, fer squeeze despres pq doni bé
            % comprovar amb el que ja hi h
        end

        function alpha = computeDamageDerivative(obj,jump,xV)
            isDerivativeZero = obj.checkIsDerivativeZero(jump);  % nDimJumpNorm (1) x ngauss x nelem
            alpha      = obj.jumpCrit .* obj.jumpFinal ./ ((obj.jumpFinal-obj.jumpCrit) .*obj.lambdaTrial.^3) ;
            % ddot       = alpha .* jump .* isDerivativeZero;
        end

        function isDerivativeZero = checkIsDerivativeZero(obj,jump) 
            obj.computeJumpNorm(jump);
            temp1 = obj.lambdaTrial - obj.jumpCrit; % f - a
            temp2 = obj.jumpFinal - obj.lambdaTrial; % b - f
            isDerivativeZero = temp1.*temp2 > 0;  % 1 x ngauss x nelem
        end       

        function B = computeMixedModeRatio(obj,jump,isJump) % jumpNorm = 1xngaussxnelem | jumpshear = 1xngaussxnelem | B = 1xngaussxnelem 
            unoZero   = ConstantFunction.create([1;0],jump.mesh);
            jumpShear = DP(jump.',unoZero);    % COMPROVAR
            B         = isJump .* (jumpShear./obj.lambdaTrial).^2;
        end

        function computeJumpLimits(obj,jump)
            isJump = obj.lambdaTrial > 1e-12;
            B = obj.computeMixedModeRatio(jump,isJump);
                B = 0;      % Impose mode 1 check
            obj.jumpCrit  = sqrt(obj.jump0Normal.^2 + (obj.jump0Shear.^2 - obj.jump0Normal.^2).*B.^obj.eta);
            obj.jumpFinal = (obj.jump0Normal*obj.jumpFinalNormal + (obj.jump0Shear*obj.jumpFinalShear - obj.jump0Normal*obj.jumpFinalNormal).*B.^obj.eta)./obj.jumpCrit;
            
            % obj.jumpCrit  = 1.25e-7; 
            % obj.jumpFinal = 0.025e-3;
        end
    end
end