classdef PhaseFieldAngleUpdater < handle
    
    properties (Access = private)
        phiOld
    end

    properties (Access = private)
        mesh
        dmgLimit
    end

    methods (Access = public)

        function obj = PhaseFieldAngleUpdater(cParams)
            obj.init(cParams);
        end

        function updatePhiOld(obj,phi)
            obj.phiOld.setFValues(phi.fValues);
        end

        function [alpha] = update(obj,u,mat,phiAlpha,alphaOld)
            s.operation = @(xV) obj.computeNewOrientation(u,mat,phiAlpha,xV);
            s.mesh = obj.mesh;
            s.ndimf = 1;
            alphaNew = DomainFunction(s);

            alpha = (obj.phiOld >= obj.dmgLimit).*alphaOld + (obj.phiOld < obj.dmgLimit).*alphaNew;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.initPhi.fun.mesh;
            obj.phiOld   = copy(cParams.initPhi.fun);
            obj.dmgLimit = 1e-20;
        end

        function alpha = computeNewOrientation(obj,u,mat,phiAlpha,xV)
            eps = SymGrad(u).evaluate(xV);         
            ex = squeeze(eps(1,1,:,:));
            ey = squeeze(eps(2,2,:,:));
            exy = squeeze(eps(1,2,:,:));

            % Perpendicular to principal directions
            % alpha = atan2(2* exy,(ex-ey))./2 + pi/2;
            % alpha = alpha';

            %Pedersen
            C = mat.obtainBaseTensor(phiAlpha);
            Cmat = C.evaluate(xV);
            C11 = squeeze(Cmat(1,1,1,1,:));
            C12 = squeeze(Cmat(1,1,2,2,:));
            C22 = squeeze(Cmat(2,2,2,2,:)); 
            C33 = squeeze(Cmat(1,2,1,2,:));

            em = (ex+ey)/2; R = sqrt(((ex-ey)/2).^2 + exy.^2);
            theta = wrapTo2Pi(atan2(2.*exy,ex-ey)./2 + (ey>ex).*pi/2);
            chi   = wrapTo2Pi(atan2(em-R,em+R) + (ey>ex).*sign(em)*pi/2 + (ex == 0 & ey == 0 & exy < 0).*pi);
            % Last condition shouldn't be necessary since tan(alpha + 180)
            % = tan(alpha). It is just written in order to properly see chi
            % state.
            deltaMat = C11 + C22 - 2.*C12 - 4.*C33;
            gamma = ((C11-C22)./deltaMat).*((1+tan(chi))./(1-tan(chi)));

            lowShear = (deltaMat > 0);
            vertexInBounds = (abs(gamma)<=1);

            
            beta = lowShear.*vertexInBounds.*(0.5*acos(-gamma)) + ...
                  (lowShear.*(~vertexInBounds)).*((gamma>1).*pi/2) + ...
                  (~lowShear).*((gamma<1).*pi/2);
            if imag(beta) < 1e-5
                beta = real(beta);
            end
            alpha = wrapTo2Pi((beta + theta)');

            if size(xV,2) == 4
            load('angleInfo','angleInfo')
            angleInfo.alpha = [angleInfo.alpha alpha'];
            angleInfo.beta  = [angleInfo.beta beta];
            angleInfo.theta = [angleInfo.theta theta];
            angleInfo.chi   = [angleInfo.chi chi];
            angleInfo.phi   = [angleInfo.phi phiAlpha.evaluate(xV)'];
            angleInfo.phiOld = [angleInfo.phiOld obj.phiOld.evaluate(xV)'];
            save('angleInfo','angleInfo')
            else
                fprintf('--------------')
                eps
                % theta
                % chi
                % vertexInBounds
                % lowShear
                % beta
                % alpha
            end
        end

    end

end