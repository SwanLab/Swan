clear all
close all
fid = fopen('executed_cases.txt');
C = textscan(fid, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
fclose(fid);
k = 1;
for i = 2:2:26
Datos{k} = [C{i}];
k = k + 1;
end
%NewData = reshape(Datos,[],5);
tipo = [Datos{1}];
Phi = Datos{2};
Theta = Datos{3};
Vol = Datos{4};
Lagrange0 = Datos{5};
Iter = Datos{6};
LagrangeInf = Datos{7};
Cost0 = Datos{8};
Cost_Inf = Datos{9};
%h_0 = Datos(:,11);
h_inf = Datos{10};
theta_0 = Datos{11};
theta_inf = Datos{12};
Energia = Datos{13};
%Residuo = Datos(:,15);


%differents_penalty = unique(Penalty);

cstring = 'rgbcmyk'; % color string


%color = {'b','r','g','m','k','y','c'};

% kpenalty = 1;
% for ipenalty = 1:length(differents_penalty)
% 
%     rows = Penalty == (differents_penalty(ipenalty)) & (abs(h_inf - Cost_Inf)./abs(h_inf) < 1*1e-6);
%     
%     if sum(rows) >= 1  
%     figure(1)    
%     h(kpenalty) = plot(Lagrange0(rows),h_inf(rows),cstring(mod(kpenalty,7)+1));
%     hold on
%     plot(Lagrange0(rows),Cost_Inf(rows),['--*',cstring(mod(kpenalty,7)+1)])
%     legend_penalty{kpenalty} = ['Penalty =  ',num2str(differents_penalty(kpenalty))];
%     kpenalty = kpenalty + 1 ;
%     legend(h,legend_penalty)
%     end
% end




% kpenalty = 1;
% for ipenalty = 1:length(differents_penalty)
% 
%     rows = Penalty == (differents_penalty(ipenalty)) & (Lagrange0 == 15) & (abs(h_inf - Cost_Inf)./abs(h_inf) < 1*1e-6);
%     
%     if sum(rows) >= 1  
%       filas = 1:length(Penalty);
%       nfila = filas(rows);
%       penaltyplot(kpenalty) = Penalty(nfila(end));
%       cost_infplot(kpenalty) = Cost_Inf(nfila(end));
%       lagrange_infplot(kpenalty) = LagrangeInf(nfila(end));
%       kpenalty = kpenalty + 1 ;
% 
%     end
% end
% figure(2)
% h2 = plot(penaltyplot,cost_infplot,'-+');  
% figure(3)
% h3 = plot(penaltyplot,lagrange_infplot,'-+');

differents_func_obj = unique(tipo);
kpenalty = 1;
for ifunc_obj = 1:length(differents_func_obj)

    rows =  abs(Phi - pi/4) <= 1e-2 & abs(Theta - 0) <= 1e-2  &  strcmp(tipo,differents_func_obj(ifunc_obj)); %& (abs(h_inf - Cost_Inf)./abs(h_inf) < 1*1e-12);
    
    
    
    %penaltyplot(kpenalty) = Penalty(nfila(end));
    cost_infplot = Cost_Inf(rows);
    lagrange_0plot = Lagrange0(rows);
    Vol_plot = Vol(rows);
    theta_plot = theta_inf(rows);
    iter_plot = Iter(rows);
    energia_plot = Energia(rows);
    
    [~,orden] = sort(Vol_plot);

    figure(1)
    hold on
    h1(ifunc_obj) = plot(Vol_plot(orden),cost_infplot(orden),['-*',cstring(mod(ifunc_obj,7)+1)]);
    title('cost function')
    figure(2)
    hold on
    h2(ifunc_obj) = plot(Vol_plot(orden),lagrange_0plot(orden),['-*',cstring(mod(ifunc_obj,7)+1)]);
    title('Lagrange')
    figure(3)
    hold on
    h3(ifunc_obj) = plot(Vol_plot(orden),theta_plot(orden),['-*',cstring(mod(ifunc_obj,7)+1)]);
    title('Theta')
    figure(4)
    hold on
    h4(ifunc_obj) = plot(Vol_plot(orden),iter_plot(orden),['-*',cstring(mod(ifunc_obj,7)+1)]);
    title('Iter')
    figure(5)
    hold on
    h5(ifunc_obj) = plot(Vol_plot(orden),energia_plot(orden),['-*',cstring(mod(ifunc_obj,7)+1)]);
    title('Energia')
   
end

% [~,orden] = sort(lagrange_0plot);
% 
% figure(4)
% h2 = plot(lagrange_0plot(orden),cost_infplot(orden),'-+');  
% title('cost function')
% figure(5)
% h3 = plot(lagrange_0plot(orden),Vol_plot(orden),'-+');  
% title('Volum')
% figure(6)
% h4 = plot(lagrange_0plot(orden),theta_plot(orden),'-+');  
% title('Theta')
% figure(7)
% h5 = plot(lagrange_0plot(orden),iter_plot(orden),'-+');  
% title('Iter')
% figure(8)
% h5 = plot(lagrange_0plot(orden),energia_plot(orden),'-+');  
% title('Energia')




