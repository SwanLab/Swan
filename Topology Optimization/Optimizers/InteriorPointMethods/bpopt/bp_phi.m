classdef bp_phi < handle
   properties (Access = public)
      phi
   end
   properties (Access = private)
      ph
      bp
      x
      xL
      xU
      s
      bL
      bU
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
          obj.bp = cParams.bp;
          obj.x = cParams.x;
          obj.xL = cParams.xL;
          obj.xU = cParams.xU;
          obj.s = cParams.s;
          obj.bL = cParams.bL;
          obj.bU = cParams.bU;
      end

      function computeObjective(obj)
         u.bp = obj.bp;
         u.x = obj.x;
         phiBase = bp_obj(u);
         phiBase.compute();
         obj.ph = phiBase.objectiveFunc;
      end

      function computePhi(obj)
         n = size(obj.x,2);
         m = size(obj.bL,2);
         for i = 1:n
            obj.ph = obj.ph - cParams.mu * (log(obj.x(i) - obj.xL(i)) + log(obj.xU(i) - obj.x(i)));
         end
         j = 0;
         for i = 1:m
            if(obj.bU(i) > obj.bL(i))
               j = j + 1;
               obj.ph = obj.ph - cParams.mu * (log(obj.s(j) - obj.bL(i)) + log(obj.bU(i) - obj.s(j)));
            end
         end
         obj.phi = obj.ph;
      end
   end
end