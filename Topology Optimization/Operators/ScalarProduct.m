classdef ScalarProduct < handle

 methods (Access = public, Static)

     function sp = computeL2(m,f,g)
         s.mesh         = m;
         s.quadType     = 'QUADRATIC';
         s.type         = 'ScalarProduct';
         int            = Integrator.create(s);
         sp = int.compute(f,g);
     end

     function sp = computeH1(m,f,g,eps)
         quadOrder = 'QUADRATIC';
         q = Quadrature.set(m.type);
         q.computeQuadrature(quadOrder);
         Df  = f.computeGradient(q);
         Dg  = g.computeGradient(q);         
         s.mesh         = m;
         s.quadType     = quadOrder;
         s.type         = 'ScalarProduct';
         int            = Integrator.create(s);
         spM = int.compute(f,g);
         spK = int.compute(Df,Dg);
         sp  = spM + eps^2*spK;
     end
   
 end
    
end
