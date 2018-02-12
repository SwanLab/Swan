function Compare_micro_gauss_point_vs_components

Case = 'BendingBeam';%'Airfoil';%'BendingBeam';

Case_A = [Case,'/',Case,'_NoMacro_Micro_Gauss'];
Case_B = [Case,'/',Case,'_NoMacro_Micro_Group'];

Path = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/';


Case_1 = Case_A;
Case_2 = Case_B;


Name =  {Case_1,Case_2};
data_plot.leyenda = {'No micro change','Micro by gauss point','Micro by components'};
data_plot.colors = 'brgcmyk';
data_plot.marker = {'*','*'};
line_style1 = {'-','--'};
line_style2 = {'-','-'};
tipo_de_iteracion = {'Ch','sigma'};
tipo_selec_micro = {'Gauss','Group'};
data_plot.linewidth = 4;


for icases = 1:length(Name)
Data{icases} = Obtaindata([Path,Name{icases}]);    
data_plot.function = Data{icases}.Cost(Data{icases}.Type == 2);
iterations(icases) = length(data_plot.function);
max_case(icases) = max(data_plot.function);
min_case(icases) = min(data_plot.function);
end



fig = randi([0 1000],1);
figure1 = figure(fig);
figure2 = figure(fig+1);

h = zeros(max(iterations),length(Name)+1);
h2 = zeros(max(iterations),length(Name)+1);


figure(figure1)
h(1,1) = plot([0 max(iterations)],[1 1],'k','LineWidth',data_plot.linewidth);

set(gca, 'FontName', 'Arial')
set(gca,'FontSize', 18)
xlabel('Iterations')
ylabel('Compliance')

%xlhand = get(gca,'xlabel');
%set(xlhand,'string','X','fontsize',20)

axis([0 1.05*max(iterations) 0.95*min(min_case) 1.05*max(max_case)])
print(figure1,[Path,Case,'/initial_compliance'],'-dpng')
hold on

figure(figure2)
hold on
h2(1,1) = plot([0 max(iterations)],[1 1],'k','LineWidth',data_plot.linewidth);

set(gca, 'FontName', 'Arial')
set(gca,'FontSize', 18)
xlabel('Iterations')
ylabel('Compliance')

axis([0 1.05*max(iterations)/2 0.95*min(min_case) 1.05*max(max_case)])
print(figure1,[Path,Case,'/initial_compliance_only_equil'],'-dpng')
hold on

for icases = 1:length(Name)
%line_style_plot = line_style{icases};
data_plot.color_plot = data_plot.colors(icases);
data_plot.function = Data{icases}.Cost(Data{icases}.Type == 2);


interval = [1:length(data_plot.function)];
h(1:length(interval)-1,icases+1) = plot_case(icases,Path,Case,interval,data_plot,tipo_selec_micro{icases},figure1,line_style1);

figure(figure2)
hold on
interval = [1:2:length(data_plot.function)];
h2(1:length(interval)-1,icases+1) = plot_case(icases,Path,Case,interval,data_plot,['only_equil',tipo_selec_micro{icases}],figure2,line_style2);


end

legend(h(1,:),data_plot.leyenda,'FontSize',10,'FontWeight','bold','Location','best')
print(figure1,[Path,Case,'/F_obj'],'-dpng')

legend(h2(1,:),data_plot.leyenda,'FontSize',10,'FontWeight','bold','Location','best')
print(figure2,[Path,Case,'/F_obj_only_equil'],'-dpng')


 % color string

end

function Data = Obtaindata(Path)
unix(['cd ',Path]);
fid = fopen([Path,'/f_obj.txt']);
tline = fgets(fid);
C = textscan(fid, '%f %f');
fclose(fid);
k = 1;
%NewData = reshape(Datos,[],5);
Data.Type = C{1};
Data.Cost = C{2};

end


function h = plot_case(icases,Path,Case,interval,data_plot,name_plot,figure_number,line_style)

figure(figure_number)
plot(0,data_plot.function(1),'LineWidth',data_plot.linewidth,'Color',data_plot.color_plot,'Marker',data_plot.marker{icases});  
print(figure_number,[Path,Case,'/',name_plot,'_',num2str(1)],'-dpng')

for iplot = 1:length(interval)-1
intervalo = [iplot iplot+1];
h(iplot,1) = plot(intervalo-1,data_plot.function(interval(intervalo)));  
line_style_plot = line_style{mod(iplot+1,2)+1};
set(h(iplot,1),'LineWidth',data_plot.linewidth,'Color',data_plot.color_plot,'Marker',data_plot.marker{icases},'LineStyle',line_style_plot)
print(figure_number,[Path,Case,'/',name_plot,'_',num2str(iplot+1)],'-dpng')
end


end



