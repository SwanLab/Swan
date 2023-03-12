classdef bp_accept < handle
    properties (Access = public)
        accept
    end
    properties (Access = private)
        gamma_phi
        gamma_th
        better_on
        filter_on
        n
        filter
        bp
        bL
        bU
        xa
        xL
        x
        s
    end

    methods (Access = public)
        function obj = bp_accept(cParams)
            obj.init(cParams);
        end

        function compute(obj)
            obj.loadParams();
            obj.checkViolations();
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.x = cParams.x;
            obj.xa = cParams.xa;
            obj.xL = cParams.xL;
            obj.xU = cParams.xU;
            obj.filter = cParams.filter;
            obj.bU = cParams.bU;
            obj.bL = cParams.bL;
            obj.s = cParams.s;
        end

        function loadParams(obj)
            obj.filter_on = false;
            obj.better_on = false;
            obj.n = size(obj.x,2);
            obj.gamma_th = 10^-5;
            obj.gamma_phi = 10^-5;
        end

        function checkViolations(obj)
            viol = 0;
            for i = 1:obj.n,
                if(obj.xa(i) < obj.xL(i)),
                    viol = viol + 1;
                end
                if(obj.xa(i) > obj.xU(i)),
                    viol = viol + 1;
                end
            end
            s.bp = obj.bp;
            s.x = obj.xa;
            s.xL = obj.xL;
            s.xU = obj.xU;
            s.bU = obj.bU;
            s.bL = obj.bL;
            s.s = obj.s;
            phiA = obj.computePhi(s);

            s.x = obj.x;
            phiX = obj.computePhi(s);

            thetaX = obj.computeTheta(s);

            s.x = s.xa;
            thetaA = obj.computeTheta(s);
            if (viol==0)
            % either condition better compared to last iteration
            if (obj.better_on)
                better = (phiA <= phiX - obj.gamma_phi*thetaX) || (thetaA <= (1 - obj.gamma_th)* thetaX);
            else
                better = true;
            end
                
            % apply filter to determine whether to accept
            if (better)
                ac = true;
                if (obj.filter_on)
                    nf = size(obj.filter,1);
                    for i = 1:nf
                        if (((1 - obj.gamma_th)*thetaA > obj.filter(i,1)) || (phiA - obj.gamma_phi*thetaA > obj.filter(i,2)))
                            ac = false;
                        end
                    end
                end
            else
                ac = false;
            end
            else
                ac = false;
            end
            obj.accept = ac;
        end
    end
    methods (Static, Access = private)
        function phi = computePhi(cParams)
            cPhi = bp_phi(cParams);
            cPhi.compute();
            phi = cPhi.phi;
        end

        function theta = computeTheta(cParams)
            cTheta = bp_theta(cParams);
            cTheta.compute();
            theta = cTheta.theta;
        end
    end
end