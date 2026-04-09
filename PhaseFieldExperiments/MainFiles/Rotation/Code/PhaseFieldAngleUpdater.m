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

        function [theta] = update(obj,u,thetaOld)
            s.operation = @(xV) obj.computeNewOrientation(u,xV);
            s.mesh = obj.mesh;
            s.ndimf = 1;
            thetaNew = DomainFunction(s);

            theta = (obj.phiOld >= obj.dmgLimit).*thetaOld + (obj.phiOld < obj.dmgLimit).*thetaNew;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.mesh     = cParams.initPhi.fun.mesh;
            obj.phiOld   = copy(cParams.initPhi.fun);
            obj.dmgLimit = 1e-12;
        end

        function theta = computeNewOrientation(~,u,xV)
            eps = SymGrad(u).evaluate(xV);         
            ex = eps(1,1,:,:);
            ey = eps(2,2,:,:);
            exy = eps(1,2,:,:);
            theta = atan2(exy,(ex-ey))./2 + pi/2;
            %(ex>ey)*pi/2 this sould be the principal directions
        end

    end

end