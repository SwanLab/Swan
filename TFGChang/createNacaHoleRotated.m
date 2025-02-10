function fV = createNacaHole(x,y,s)
            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = -deg2rad(s.AoA);
        
            x0    = s.xLE;
            y0    = s.yLE/c;
            xNaca = (x-x0)/c;
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>=p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1015*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>=p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = y0 + yc + yt.*cos(theta);
            yl    = y0 + yc - yt.*cos(theta);
        
            xu = (xNaca>=0 & xNaca<=1).*(xNaca).*cos(AoA) - (yu - y0).*sin(AoA) + x0/c;
            yu = (xNaca>=0 & xNaca<=1).*(xNaca).*sin(AoA) + (yu - y0).*cos(AoA) + y0;
            xl = (xNaca>=0 & xNaca<=1).*(xNaca).*cos(AoA) - (yl - y0).*sin(AoA) + x0/c;
            yl = (xNaca>=0 & xNaca<=1).*(xNaca).*sin(AoA) + (yl - y0).*cos(AoA) + y0;
        
            x0Rotated    = min([xu(:); xl(:)]);
            xNacaRotated = (x-x0Rotated)/(c*cos(AoA));
            f(:,:,:,1)   = yl*c - y;
            f(:,:,:,2)   = y - yu*c;
            f(:,:,:,3)   = xNacaRotated - 1; 
            f(:,:,:,4)   = -xNacaRotated;
            fV           = -max(f,[],4);        
end