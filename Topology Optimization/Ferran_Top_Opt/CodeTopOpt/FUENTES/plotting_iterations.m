function plotting_iterations(vcost)
filas = (vcost(:,1) ~= 0);
nfilas = length(vcost(filas,1));
iterations = 0;
v_opt = [];
it_plot = [];
for ifilas = 2:nfilas
 columnas = (vcost(ifilas,:) ~= 0);
 ncolumnas = length(vcost(ifilas,columnas));
 iterations = [iterations(end):iterations(end)+ncolumnas-1];
 v_plot = vcost(ifilas,columnas);
 v_opt = [v_opt v_plot(end)];
 it_plot = [it_plot iterations(end)];
 figure(10)
 hold on
 plot(iterations,v_plot,'b')
 hold on
end
figure(10)
hold on
plot(it_plot,v_opt,'r')
hold off
end