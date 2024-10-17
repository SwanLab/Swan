clear
close all

M=9/100;
p=4/10;
t=12/100;
pas=0.001;
x_p=[0:pas:1]; %S'ha de retallar una mica la punta perquè sinó queden els munts malament cap al caire de sortida 

yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1015*x_p.^4);

for j=1:1:size(x_p,2)
    if x_p(j)<=p
        y_c(j)=(M/(p^2))*(2*p*x_p(j)-x_p(j)^2);
    elseif x_p(j)>p
        y_c(j)=(M/(1-p)^2)*((1-2*p)+2*p*x_p(j)-x_p(j)^2);
    end
end

%Plot chamber line:
plot(x_p,y_c)
axis equal
grid on
hold on
% scatter(x_p,y_c);
hold on
%Plot airfoil with circles:
for ii=1:1:size(x_p,2)
    x_c = [x_p(ii)-yt(ii):0.0001:x_p(ii)+yt(ii)+0.0001];
    y = sqrt(yt(ii)^2 - (x_c-x_p(ii)).^2);


    plot(x_c,y+y_c(ii),'b');
    hold on
    plot(x_c,-y+y_c(ii),'b');
    hold on

end

axis equal

% %% Nou mètode
% M=0/100;
% p=4/10;
% t=12/100;
% 
% pas=0.00001;
% 
% x_p=[0:pas:1];
% 
% yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1015*x_p.^4); % el valor que hi ha ara és per tenir la punta tancada-0.1036
% 
% for j=1:1:size(x_p,2)
%     if x_p(j)<=p
%         y_c(j)=(M/(p^2))*(2*p*x_p(j)-x_p(j)^2);
%         theta(j)=atan((2*M/(p^2))*(p-x_p(j)));
%     elseif x_p(j)>p
%         y_c(j)=(M/(1-p)^2)*((1-2*p)+2*p*x_p(j)-x_p(j)^2);
%         theta(j)=atan(((2*M/((1-p)^2))*(p-x_p(j))));
%     end
% end
% 
% x_u = x_p - yt.*sin(theta);
% y_u = y_c + yt.*cos(theta);
% x_l = x_p + yt.*sin(theta);
% y_l = y_c - yt.*cos(theta);
% 
% % for i=2:1:lenght(x_u)
% %Plot chamber line:
% plot(x_p,y_c)
% axis equal
% grid on
% hold on
% scatter(x_p,y_c);
% hold on
% %Plot airfoil with circles:
% for ii=1:1:size(x_p,2)
%     x_c = [x_p(ii)-yt(ii):0.0001:x_p(ii)+yt(ii)+0.0001];
%     y = sqrt(yt(ii)^2 - (x_c-x_p(ii)).^2);
% 
% 
%     plot(x_c,y+y_c(ii));
%     hold on
%     plot(x_c,-y+y_c(ii));
%     hold on
% 
% end
% 
% axis equal

%% Mètode de les circumf
clear

M=9/100;
p=4/10;
t=12/100;

pas=0.00001;

x_p=[0:pas:1];

yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1015*x_p.^4); % el valor que hi ha ara és per tenir la punta tancada-0.1036

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

plot(x_coord,y_coord,'r','LineWidth',2);
hold on
plot(x_l,y_l);

axis equal

% %% GPT
% clear
% 
% % Paràmetres del perfil NACA 4 dígits
% m = 5/100;    % coeficient de curvatura màxima
% p = 4/10;     % posició de la curvatura màxima (en x/c)
% t = 24/100;    % espessor màxim (en t/c)
% c = 1;       % longitud de la corda
% 
% 
% % Nombre de punts
% N = 100;
% 
% % Coordenades x
% x = linspace(0, c, N);
% 
% % Calcula el gruix a cada punt x
% yt = (t/0.2)*c*(0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c).^2 + 0.2843*(x/c).^3 - 0.1015*(x/c).^4);
% 
% % Calcula la línia mitjana de curvatura yc i la seva derivada dyc/dx
% yc = zeros(1, N);
% dyc_dx = zeros(1, N);
% 
% for i = 1:N
%     if x(i) <= p*c
%         yc(i) = (m/(p^2))*(2*p*(x(i)/c) - (x(i)/c)^2);
%         dyc_dx(i) = (2*m/(p^2))*(p - x(i)/c);
%     else
%         yc(i) = (m/((1-p)^2))*(1 - 2*p + 2*p*(x(i)/c) - (x(i)/c)^2);
%         dyc_dx(i) = (2*m/((1-p)^2))*(p - x(i)/c);
%     end
% end
% 
% % Angle theta en cada punt x
% theta = atan(dyc_dx);
% 
% % Coordenades del perfil superior i inferior
% xu = x - yt.*sin(theta);
% yu = yc + yt.*cos(theta);
% xl = x + yt.*sin(theta);
% yl = yc - yt.*cos(theta);
% 
% % Dibuixa el perfil NACA
% figure;
% hold on;
% plot(xu, yu, 'k', 'LineWidth', 1.5);
% plot(xl, yl, 'k', 'LineWidth', 1.5);
% 
% % Dibuixa els cercles al llarg de la línia mitjana de curvatura
% for i = 1:N
%     radius = yt(i) / 2; % Radi del cercle és la meitat del gruix en aquest punt
%     theta_circle = linspace(0, 2*pi, 50);
%     xc = x(i) + radius * cos(theta_circle); % Coordenades x dels cercles
%     yc_circle = yc(i) + radius * sin(theta_circle); % Coordenades y dels cercles
%     plot(xc, yc_circle, 'r'); % Ploteja els cercles
% end
% 
% axis equal;
% title('Perfil NACA 4 dígits amb cercles sobre la línia mitjana de curvatura');
% xlabel('x');
% ylabel('y');
% hold off;


%% inpolygon
% x(1,:,:)

% x=[0,0.5,0.8,2];
% y=[0,0,0.01,-2];
% 
% 
% GG =  double(inpolygon(x,y,x_coord,y_coord)); %@(x)
% 
% [in,on] =  inpolygon(x,y,x_coord,y_coord); %@(x)
















