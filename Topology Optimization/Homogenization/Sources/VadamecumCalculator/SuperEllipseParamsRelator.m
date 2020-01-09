classdef SuperEllipseParamsRelator < handle
    
    
   methods (Access = public, Static)
       
       function mx = mx(txi,rho,q)
            c = SuperEllipseParamsRelator.c();
            n = (1-rho)*tan(txi);
            d = c(q);
            mx = sqrt(n/d);               
       end
       
       function my = my(txi,rho,q)
            c = SuperEllipseParamsRelator.c();
            n = (1-rho);
            d = c(q)*tan(txi);
            my = sqrt(n/d);                      
       end
       
       function rho = rho(mx,my,q)
            c = SuperEllipseParamsRelator.c();
            rho =  1 - c(q)*mx*my;
       end
       
       function txi = txi(mx,my)
           txi = atan(mx/my);
       end
       
       function c = c()       
          c = @(q) gamma(1 + 1/q)^2/gamma(1 + 2/q); 
       end
       
       function my = myFromMxAndTxi(mx,txi)
          my = mx./tan(txi);
       end
       
       function mx = mxFromMyAndTxi(my,txi)
          mx = my.*tan(txi);
       end
       
       function rho = rhoFromMxAndTxi(mx,txi,q)
           sE  = SuperEllipseParamsRelator();
           my  = sE.myFromMxAndTxi(mx,txi);
           rho = sE.rho(mx,my,q);
       end
       
       function rho = rhoFromMyAndTxi(my,txi,q)
           sE  = SuperEllipseParamsRelator();
           mx  = sE.mxFromMyAndTxi(my,txi);
           rho = sE.rho(mx,my,q);
       end       
       
       function txi = txiFromMxAndRho(mx,rho,q)
           sE  = SuperEllipseParamsRelator();
           c  = sE.c();
           txi = atan(c(q)*mx^2/(1-rho));
       end
       
       function txi = txiFromMyAndRho(my,rho,q)
           sE  = SuperEllipseParamsRelator();
           c  = sE.c();
           txi = atan((1-rho)/(c(q)*my^2));
       end       
       
   end
    
    
end