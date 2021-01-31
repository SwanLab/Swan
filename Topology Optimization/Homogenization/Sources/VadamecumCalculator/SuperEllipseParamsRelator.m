classdef SuperEllipseParamsRelator < handle
       
   methods (Access = public, Static)
       
       function mx = mx(xi,rho,q)
            c = SuperEllipseParamsRelator.c();
            n = (1-rho).*tan(xi);
            d = c(q);
            mx = sqrt(n./d);               
       end
       
       function my = my(xi,rho,q)
            c = SuperEllipseParamsRelator.c();
            n = (1-rho);
            d = c(q).*tan(xi);
            my = sqrt(n./d);                      
       end
       
       function rho = rho(mx,my,q)
            c = SuperEllipseParamsRelator.c();
            rho =  1 - c(q).*mx.*my;
       end
       
       function xi = xi(mx,my)
           xi = atan(mx./my);
       end
       
       function c = c()       
          c = @(q) gamma(1 + 1./q).^2./gamma(1 + 2./q); 
       end
       
       function my = myFromMxAndTxi(mx,xi)
          my = mx./tan(xi);
       end
       
       function mx = mxFromMyAndTxi(my,xi)
          mx = my.*tan(xi);
       end
       
       function rho = rhoFromMxAndTxi(mx,xi,q)
           sE  = SuperEllipseParamsRelator();
           my  = sE.myFromMxAndTxi(mx,xi);
           rho = sE.rho(mx,my,q);
       end
       
       function rho = rhoFromMyAndTxi(my,xi,q)
           sE  = SuperEllipseParamsRelator();
           mx  = sE.mxFromMyAndTxi(my,xi);
           rho = sE.rho(mx,my,q);
       end       
       
       function xi = xiFromMxAndRho(mx,rho,q)
           sE  = SuperEllipseParamsRelator();
           c  = sE.c();
           xi = atan(c(q).*mx.^2./(1-rho));
       end
       
       function xi = xiFromMyAndRho(my,rho,q)
           sE  = SuperEllipseParamsRelator();
           c  = sE.c();
           xi = atan((1-rho)./(c(q).*my.^2));
       end 
       
       function xiUB = xiUB(rho)
           sE  = SuperEllipseParamsRelator();   
           qMax = 32;
           qMin = 2;
           mxUB = 0.99;
           myLB = 0.01;
           xi1 = sE.xiFromMxAndRho(mxUB,rho,qMax);
           xi2 = sE.xiFromMyAndRho(myLB,rho,qMin);           
           xiUB = min(xi1,xi2);
       end

       function xiLB = xiLB(rho)
           sE  = SuperEllipseParamsRelator();   
           qMax = 32;
           qMin = 2;
           mxLB = 0.01;                      
           myUB = 0.99;           
           xi1 = sE.xiFromMxAndRho(mxLB,rho,qMin);
           xi2 = sE.xiFromMyAndRho(myUB,rho,qMax);           
           xiLB = min(xi1,xi2);
       end
       
       function [mXb,mYb] = updateMxMyWithOtherQ(mXa,mYa,q,qNew)
            sE   = SuperEllipseParamsRelator;
            rho  = sE.rho(mXa,mYa,q);
            xi   = sE.xi(mXa,mYa);           
            mXb  = sE.mx(xi,rho,qNew);
            mYb  = sE.my(xi,rho,qNew);
       end
       
       
   end
    
    
end