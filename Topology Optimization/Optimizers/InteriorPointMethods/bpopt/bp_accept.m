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
        end

        function loadParams(obj)
            obj.filter_on = false;
            obj.better_on = false;
            obj.n = size(obj.x,2);
            obj.gamma_th = 10^-5;
            obj.gamma_phi = 10^-5;
        end

        function phi = computePhi(cParams)  % pasar a estatica

        end

        function theta = computeTheta(obj)          % pasar a estática
        
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
            s.xa = obj.xa;
            s.xL = obj.xL;
            s.xU = obj.xU;
            phiA = obj.computePhi(s);

            c.bp = obj.bp;
            c.x = obj.x;
            c.xL = obj.xL;
            c.xU = obj.xU;
            phiX = obj.computePhi(c);

            w.bp = obj.bp;
            w.x = obj.x;
            thetaX = obj.computeTheta(w);

            u.bp = obj.np;
            u.xa = obj.xa;
            thetaA = obj.computeTheta(u);
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
end

function [ac] = bp_accept(bp,x,xa,xL,xU,filter)
% parameters
filter_on = false;
better_on = false;
n = size(x,2);
gamma_th = 10^-5;
gamma_phi = 10^-5;

% check for violations of variable constraints
% shouldn't need to check these
viol = 0;
for i = 1:n,
    if(xa(i) < xL(i)),
        viol = viol + 1;
    end
    if(xa(i)> xU(i)),
        viol = viol + 1;
    end
end

if (viol==0),
  % either condition better compared to last iteration
  if (better_on),
     better = (bp_phi(bp,xa,xL,xU) <= bp_phi(bp,x,xL,xU) - gamma_phi*bp_theta(bp,x)) || (bp_theta(bp,xa) <= (1-gamma_th)* bp_theta(bp,x));
  else
     better = true;
  end
     
  % apply filter to determine whether to accept
  if (better),
      ac = true;
      if (filter_on),
         nf = size(filter,1);
         for i = 1:nf,
             if (((1 - gamma_th)*bp_theta(bp,xa) > filter(i,1)) || (bp_phi(bp,xa,xL,xU) - gamma_phi*bp_theta(bp,xa) > filter(i,2))),
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