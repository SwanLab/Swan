        function fV = createNacaHole(x,y,s)
            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
            x0     = s.xLE;
            y0     = s.yLE/c;
            offsetX  = (x - x0)/c;
            offsetY  = y/c - y0;

            xNaca    = offsetX.*cos(AoA) - offsetY.*sin(AoA);
            yNaca    = offsetX.*sin(AoA) + offsetY.*cos(AoA);
        
            yc   = (xNaca>=0 & xNaca<=p).*(m./p^2.*(2*p*xNaca-xNaca.^2))+...
                    (xNaca>=p & xNaca<=1).*(m./(1-p)^2.*((1-2*p)+2*p*xNaca-xNaca.^2));
            yt   = (xNaca>=0 & xNaca<=1).*(5*t*(0.2969*sqrt(xNaca)-0.1260*xNaca-0.3516*xNaca.^2+0.2843*xNaca.^3-0.1015*xNaca.^4));
            dydx = (xNaca>=0 & xNaca<=p).*(2*m/p^2.*(p-xNaca))+...
                    (xNaca>=p & xNaca<=1).*(2*m/(1-p)^2.*(p-xNaca));
        
            theta = atan(dydx);
            yu    = yc + yt.*cos(theta);
            yl    = yc - yt.*cos(theta);

            f(:,:,:,1)   = yl - yNaca;
            f(:,:,:,2)   = yNaca - yu;
            f(:,:,:,3)   = xNaca - 1; 
            f(:,:,:,4)   = -xNaca;

            % No se quina d'aquestes dues funcionaria. Jo primer resoldria
            % el tema de la rotació i després revisem els detallets dels
            % forats petits

            %f  = f./max(abs(f));

            % f(:,:,:,1) = f(:,:,:,1)./max(abs(f(:,:,:,1)),[],'all');
            % f(:,:,:,2) = f(:,:,:,2)./max(abs(f(:,:,:,2)),[],'all');
            % f(:,:,:,3) = f(:,:,:,3)./max(abs(f(:,:,:,3)),[],'all');
            % f(:,:,:,4) = f(:,:,:,4)./max(abs(f(:,:,:,4)),[],'all');

            fV = -max(f,[],4);        
        end


function fV = createNacaHole(x,y,s)
            c   = s.chord;
            p   = s.p;
            m   = s.m;
            t   = s.t;
            AoA = deg2rad(s.AoA);
        
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

            f(:,:,:,1) = yl*c - y;
            f(:,:,:,2) = y - yu*c;
            f(:,:,:,3) = xNaca - 1; 
            f(:,:,:,4) = -xNaca;
            fV         = -max(f,[],4);  
        end