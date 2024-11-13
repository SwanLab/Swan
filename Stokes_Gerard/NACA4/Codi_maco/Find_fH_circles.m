function fH = Find_fH_circles(M,p,t,x_centr,y_centr,AOAd)

pas=0.001;

x_p=[0:pas:1-pas*15]; %S'ha de retallar una mica la punta perquè sinó queden els munts malament cap al caire de sortida 

yt = 5*t*(0.2969*sqrt(x_p)-0.1260*x_p-0.3516*x_p.^2+0.2843*x_p.^3-0.1015*x_p.^4);

% for j=1:1:size(x_p,2)
%     if x_p(j)<=p
%         y_c(j)=(M/(p^2))*(2*p*x_p(j)-x_p(j)^2);
%     elseif x_p(j)>p
%         y_c(j)=(M/(1-p)^2)*((1-2*p)+2*p*x_p(j)-x_p(j)^2);
%     end
% end

y_c = (  (M/(p^2))*(2*p*x_p-x_p.^2)  ).*(x_p<=p) + ...
        (    (M/(1-p)^2)*((1-2*p)+2*p*x_p-x_p.^2)  ).*(x_p>p);

% Plot chamber line:
plot(x_p,y_c)
axis equal
grid on

%Plot airfoil with circles:
figure
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

% minxC = min(x_p-yt);
% maxxC = max(x_p+yt);
% x_c   = linspace(minxC,maxxC,length(yt));
% y     = sqrt(yt.^2 - (x_c-x_p).^2);
% plot(x_c,y+y_c);
% hold on
% plot(x_c,-y+y_c);

axis equal

x_le = x_centr-0.5;
AOA = -deg2rad(AOAd);

rn  = yt;
x_cnr = x_p+x_le;
y_cnr = y_c+y_centr;

x_cn = (x_cnr-x_centr).*cos(AOA)-(y_cnr-y_centr).*sin(AOA)+x_centr;
y_cn = (x_cnr-x_centr).*sin(AOA)+(y_cnr-y_centr).*cos(AOA)+y_centr;

terms = cell(1, length(rn));

for jj = 1:length(rn)
    terms{jj} = sprintf('((x(1,:,:)-%f).^2 + (x(2,:,:)-%f).^2 - %f.^2)', x_cn(jj), y_cn(jj), rn(jj)); %Crea els termes per cada cercle (per després ajuntar-ho)
end

while length(terms) > 1
    new_terms = {};
    for jj = 1:2:length(terms)-1
        new_terms{end+1} = sprintf('min(%s, %s)', terms{jj}, terms{jj+1}); %Anem creant parelles de termes a dins del min()
    end
    if mod(length(terms), 2) == 1
        new_terms{end+1} = terms{end}; %Si queda algun terme sol, el posem al final
    end
    terms = new_terms; %Assignem els termes i tornem a mirar si n'hi ha més d'un. Al final, tots s'hauran anat ajuntat fins només quedar 1.
end

func_str = ['@(x) -', terms{1}];
fH = str2func(func_str);






end




















