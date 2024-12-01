clear all
format compact

%% Part extra (planejador)

% Els models d'avions utilitzats són extrets d'exercicis de l'assignatura 
% 'Vehicles Aeroespacials' de la carrera. 

g=9.81;
S=12.5; b=15; Al=b^2/S; e=0.95;
cd0=0.01; 
E=@(cl) cl./(cd0+cl.^2./(pi*Al*e)); % Eficiència aerodinàmica
dE=@(cl) cd0+cl.^2./(pi*Al*e) - 2.*cl.^2/(pi*Al*e); % Eficiència màxima
[xk,res,it] = newton(0.5,1e-8,15,dE,'0'); % Newton per resoldre 
cl=xk(end) 
Ef=E(cl)
cd=cl/Ef;
lam=-0.0065; R=287; T0=288.16; p0=1.225; % Atmosfera estàndard

M=0.1; a=sqrt(R*1.4*(T0+2000*lam)); % Mach i velocitat del so

T=@(h) T0+lam.*h;
p=@(h) p0.*(T(h)./T0).^(-g./(R.*lam)-1);
D=@(v,h) 0.5.*v.^2.*p(h).*S.*cd;
L=@(v,h) 0.5.*v.^2.*p(h).*S.*cl;
Dra0=D(M*a,2000); % Drag al moment 0
Lif0=L(M*a,2000); % Lift al moment 0 
W=600*g; 
planeo=@(t,q) [q(3).*cos(q(4)); q(3).*sin(q(4)); g.*(-sin(q(4))-...
    D(q(3),q(2))./W); g./q(3).*(L(q(3),q(2))./W-cos(q(4)))]; 
% Sistema d'equacions que defineix el problema 

tspan=[0 2162]; % Interval de temps, utilitzant ja com a temps final el 
% valor que es calcula seguidament per fer les gràfiques fins que cau

condicions=[0 2000 M*a 0];
[t,q]=ode45(planeo,tspan,condicions);

AngLim=q(end,4)*180/pi; % Valor utilitzat per estudiar la convergència 
% de l'angle
AngTeo=atan(1/Ef)*180/pi; % Valor teòric a què hauria de convergir

figure(5)
plot(q(:,1),q(:,2))
title('Trajectory of the airplane')
xlabel('Distance(m)')
ylabel('Height(m)')
axis([0 q(end,1) 0 2000])

figure(6)
plot(t/60,3.6*q(:,3))
title('Evolution of speed with time')
xlabel('Time (minutes)')
ylabel('Speed (km/h)')

figure(7)
plot(t/60,180/pi*q(:,4))
title('Evolution of angle with time')
xlabel('Time (minutes)')
ylabel('Angle (degrees)')

figure(8)
plot(180/pi*q(:,4),q(:,3))
title('Evolution of speed with angle')
xlabel('Angle (degrees)')
ylabel('Speed (km/h)')



tend=2150;
q(end,2)=1;
tspan=[0 tend];
while q(end,2)>0
    tspan=[0 tend];
    condicions=[0 2000 M*a 0];
[t,q]=ode45(planeo,tspan,condicions);
tend=tend+1;
end
time=tend % Temps fins caure al terra 
dist=q(end,1) % Distància recorreguda
vel=q(end,3)*3.6 %Velocitat en caure (km/h)


%% Part extra (avió comercial)

g=9.81;
S=428; b=61; Al=b^2/S; e=0.8;
cd0=0.02;
E=@(cl) cl./(cd0+cl.^2./(pi*Al*e)); % Eficiència
dE=@(cl) cd0+cl.^2./(pi*Al*e) - 2.*cl.^2/(pi*Al*e); % Eficiència màxima
[xk,res,it] = newton(0.5,1e-8,15,dE,'0'); 
cl=xk(end)
Ef=E(cl)
cd=cl/Ef;

M=0.6; a=sqrt(R*1.4*(T0+2000*lam)); % Mach i v_so
lam=-0.0065; R=287; T0=288.16; p0=1.225; % Atm estàndard
T=@(h) T0+lam.*h;
p=@(h) p0.*(T(h)./T0).^(-g./(R.*lam)-1);
D=@(v,h) 0.5.*v.^2.*p(h).*S.*cd;
L=@(v,h) 0.5.*v.^2.*p(h).*S.*cl;
Li0=L(M*a,2000);
Dr0=D(M*a,2000);
OEW=245000*g; PL=60000*g; FW0=125000*g;
W=OEW+PL+0.98*FW0;
planeo=@(t,q) [q(3).*cos(q(4)); q(3).*sin(q(4)); g.*(-sin(q(4))-...
    D(q(3),q(2))./W); g./q(3).*(L(q(3),q(2))./W-cos(q(4)))];
tspan=[0 271]; % Temps final agafat del resultat posterior
condicions=[0 2000 M*a pi/60];
[ta,p]=ode45(planeo,tspan,condicions);

figure(9)
plot(q(:,1),q(:,2)); hold on
plot(p(:,1),p(:,2))
title('Trajectory of airplane next to glider')
legend( 'Glider trajectory','Airplane trajectory')
xlabel('Distance(m)')
ylabel('Height(m)')

figure(10)
plot(t/60,3.6*q(:,3)); hold on
plot(ta/60,3.6*p(:,3))
title('Evolution of speed with time')
legend('Glider speed','Airplane speed')
xlabel('Time (minutes)')
ylabel('Speed (km/h)')

figure(11)
plot(q(:,1),180/pi*q(:,4)); hold on
plot(p(:,1),180/pi*p(:,4))
title('Evolution of angle with distance')
legend('Glider angle','Airplane angle')
xlabel('Distance (m)')
ylabel('Angle (degrees)')



tend=200;
p(end,2)=1;
tspan=[0 tend];
while p(end,2)>0
    tspan=[0 tend];
    condicions=[0 2000 M*a 0];
[t,p]=ode45(planeo,tspan,condicions);
tend=tend+1;
end
time=tend % Temps que tarda a caure
dist=p(end,1) % Distància recorreguda
vel=p(end,3)*3.6 % Velocitat quan cau
LandAngle=p(end,4)*180/pi % Angle amb el terra quan cau

%% Funcions

function [xk,res,it] = newton(a,tol,itmax,fun,dfun)
xk = [a]; fk = fun(xk); res = abs(fk); it = 0;
tolk = res(1); dx = 1e-8 ;
while it < itmax && tolk > tol
if dfun == '0'
dfk = (fun(xk(end)+dx)-fk)/dx;
else
dfk = feval(dfun,xk(end));
end
xk = [xk, xk(end) - fk/dfk]; fk = fun(xk(end));
res = [res abs(fk)]; tolk = abs(xk(end)-xk(end-1));
it = it + 1;
end
end