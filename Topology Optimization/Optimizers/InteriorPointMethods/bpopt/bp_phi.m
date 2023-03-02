classdef bp_phi < handle
   properties (Access = public)
      phi
   end
   properties (Access = private)
      ph
   end

   methods (Access = public)
      function obj = bp_phi(cParams)
         obj.init(cParams);
      end
      function compute(obj)
         obj.computeObjective();
         obj.computePhi();
      end
   end

   methods (Access = private)
      function init(obj,cParams)
      end

      function computeObjective(obj)
         phiBase = bp_obj(obj);
         phiBase.compute();
         obj.ph = phiBase.objectiveFunc;
      end

      function computePhi(obj)
         n = size(x,2);
         m = size(bL,2);
         for i = 1:n,
            obj.ph = obj.ph - cParams.mu * (log(x(i)-xL(i)) + log(xU(i)-x(i)));
         end
         j = 0;
         for i = 1:m,
            if(bU(i)>bL(i)),
               j = j + 1;
               obj.ph = obj.ph - cParams.mu * (log(s(j)-bL(i)) + log(bU(i)-s(j)));
            end
         end
         obj.phi = obj.ph;
      end
   end
end
function [ph] = bp_phi(bp,x,xL,xU,s,bL,bU)
    ph = bp_obj(bp,x);
    n = size(x,2);
    m = size(bL,2);
    for i = 1:n,
       ph = ph - bp.mu * (log(x(i)-xL(i)) + log(xU(i)-x(i)));
    end
    j = 0;
    for i = 1:m,
       if(bU(i)>bL(i)),
          j = j + 1;
          ph = ph - bp.mu * (log(s(j)-bL(i)) + log(bU(i)-s(j)));
       end
    end
