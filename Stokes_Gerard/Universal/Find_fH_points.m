function fH = Find_fH_points(M,p,t,x_centr,y_centr,AOAd)

pas=0.0001;

x_p=[0:pas:0.95];

yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1036*x_p.^4); %-0.1015 el valor que hi ha ara és per tenir la punta tancada

for j=1:1:size(x_p,2)
    if x_p(j)<=p
        y_c(j)=(M/(p^2))*(2*p*x_p(j)-x_p(j)^2);
        theta(j)=atan((2*M/(p^2))*(p-x_p(j)));
    elseif x_p(j)>p
        y_c(j)=(M/(1-p)^2)*((1-2*p)+2*p*x_p(j)-x_p(j)^2);
        theta(j)=atan(((2*M/((1-p)^2))*(p-x_p(j))));
    end
end

x_u = x_p - yt.*sin(theta);
y_u = y_c + yt.*cos(theta);
x_l = x_p + yt.*sin(theta);
y_l = y_c - yt.*cos(theta);

x_coord = [x_u flip(x_l)];
y_coord = [y_u flip(y_l)];

plot(x_coord,y_coord);
hold on
plot(x_l,y_l);
axis equal


x_LE=x_centr-0.5;
AOA = -deg2rad(AOAd);
x_cnr = x_coord+x_LE;
y_cnr = y_coord+y_centr;

x_cn = (x_cnr-x_centr).*cos(AOA)-(y_cnr-y_centr).*sin(AOA)+x_centr;
y_cn = (x_cnr-x_centr).*sin(AOA)+(y_cnr-y_centr).*cos(AOA)+y_centr;

fH = @(x)  double(inpolygon(x(1,:,:),x(2,:,:),x_cn,y_cn))-0.5; 





end






