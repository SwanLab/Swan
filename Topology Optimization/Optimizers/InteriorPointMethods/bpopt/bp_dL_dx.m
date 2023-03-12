classdef bp_dL_dx < handle
   properties (Access = public)
      dLagrangian
   end
   properties (Access = private)
      LBase
      L1
      sp
      xp
      x
      bL
      bU
   end

   methods (Access = public)
      function obj = bp_dL_dx(cParams)
         obj.init(cParams);
      end
      function compute(obj)
         obj.computeLagrangianBase();
         obj.computeLagrangianDerivative();
      end
   end
   methods (Access = private)
      function init(obj,cParams)
         obj.x = cParams.x;
         obj.xL = cParams.xL;
         obj.s = cParams.s;
         obj.lam = cParams.lam;
         obj.bL = cParams.bL;
         obj.bU = cParams.bU;
      end

      function computeLagrangianBase(obj)
         LB = bp_lagrangian(obj);
         LB.compute();
         obj.LBase = LB.Lagrangian;
      end

      function computeLagrangian(obj)
         L = bp_lagrangian(obj);
         L.compute();
         obj.L1 = L.Lagrangian;
      end

      function computeLagrangianDerivative(obj)
         n = size(obj.x,2);
         m = size(obj.bL,2);
         ep = 1e-5;
         for i = 1:n
            obj.xp = obj.x;
            obj.xp(i) = obj.x(i) + ep;
            obj.computeLagrangian();
            dL(i) = [(obj.L1 - obj.LBase)/ep];
         end
         k = 0;
         for i = 1:m
            if (obj.bU(i) > obj.bL(i))
               k = k + 1;
               obj.sp = onj.s;
               obj.sp(i) = s(i) + ep;
               obj.computeLagrangian();
               dL(k + n) = [(obj.L1 - obj.LBase)/ep];
            end
         end
         obj.dLagrangian = dL;
      end
   end
end