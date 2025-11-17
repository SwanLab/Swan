classdef SimpInterpolationP3Anisotropic < handle
    
   properties (Access = private)
        C1
        C0
        pExp
   end

    methods (Access = public)
        function obj = SimpInterpolationP3Anisotropic(cParams)
            obj.init(cParams)
        end

        function C = computeConsitutiveTensor(obj,rho,xV)
            m = obj.createMaterial(rho{1});
            C = m.evaluate(xV,obj.C1);
            rhoEv = rho{1}.evaluate(xV);
            nGauss = size(rhoEv,2);
            nElem = size(rhoEv,3);
            rhoEv = reshape(rhoEv,[1 1 1 1 nGauss nElem]);
            C = (1-rhoEv.^3)*1e-3.*C + (rhoEv.^3).*C;
        end

        function dC = computeConsitutiveTensorDerivative(obj,rho,xV)
            m = obj.createMaterial(rho{1});
            C = m.evaluate(xV,obj.C1);
            rhoEv = rho{1}.evaluate(xV);
            nGauss = size(rhoEv,2);
            nElem = size(rhoEv,3);
            rhoEv = reshape(rhoEv,[1 1 1 1 nGauss nElem]);
            dC = (-3.*rhoEv.^2)*1e-3.*C + 3.*(rhoEv.^2).*C; 
        end


           
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.C1 = cParams.C1;
            obj.C0 = cParams.C0;
            obj.pExp = 3;
        end

         function m = createMaterial(obj,rho)
            s.type    = 'ANISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.mesh    = rho.mesh;
            m = Material.create(s);
        end
    end
end