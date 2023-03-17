classdef LagrangianDerivativeComputer < handle
   properties (Access = public)
      dLagrangian
   end
   properties (Access = private)
      LBase
      L1
      sp
      xp
      x
      s
      bp
      bL
      bU
      lambda
   end

   methods (Access = public)
      function obj = LagrangianDerivativeComputer(cParams)
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
         obj.bp = cParams.bp;
         obj.s = cParams.s;
         obj.lambda = cParams.lam;
         obj.bL = cParams.bL;
         obj.bU = cParams.bU;
      end

      function computeLagrangianBase(obj)
         u.bp = obj.bp;
         u.x = obj.x;
         u.s = obj.s;
         u.lam = obj.lambda;
         u.bL = obj.bL;
         u.bU = obj.bU;
         LB = LagrangianComputer(u);
         LB.compute();
         obj.LBase = LB.lagrangian;
      end

      function computeLagrangian(obj)
         u.bp = obj.bp;
         u.x = obj.x;
         u.s = obj.s;
         u.lam = obj.lambda;
         u.bL = obj.bL;
         u.bU = obj.bU;
         L = LagrangianComputer(u);
         L.compute();
         obj.L1 = L.lagrangian;
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
               obj.sp = obj.s;
               obj.sp(i) = obj.s(i) + ep;
               obj.computeLagrangian();
               dL(k + n) = [(obj.L1 - obj.LBase)/ep];
            end
         end
         obj.dLagrangian = dL;
      end
   end
end