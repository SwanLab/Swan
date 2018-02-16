function CompareMacroChange_micro_gauss_point_vs_components

cost_classical_industrial_abs = 1.906210178151383e+02;
cost_total_material_abs = 57.896666138437929;
cost_classical = cost_classical_industrial_abs/cost_total_material_abs;


Case = 'Rib';%'CantileverBeam';%'Airfoil';%'BendingBeam';
Name = 'Rib_MicroRareShear_to0_1eta0_999'; %'CantileverBeam'

Case_A = [Case,'/',Name,'MacroOnly'];
Case_B = [Case,'/',Name,'MacroMicroGauss'];
Case_C = [Case,'/',Name,'MacroMicroGroup'];

Path = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/';


Case_1 = Case_A;
Case_2 = Case_B;
Case_3 = Case_C;


Name =  {Case_1,Case_2,Case_3};
data_plot.leyenda = {'No Micro Design','Continuum Micro Design','Discrete Micro Design'};
data_plot.leyenda2 = {'Continuum','Discrete'};
data_plot.colors = 'kbrgcmyk';
data_plot.marker = {'*','*','*'};
line_style1 = {'-','-','-'};
tipo_selec_micro = {'NoMicroDesign','ContinuumMicroDesign','DiscreteMicroDesign'};
data_plot.linewidth = 4;


for icases = 1:length(Name)
Data{icases} = Obtaindata([Path,Name{icases}]);    
data_plot.function = Data{icases}.Cost;
iterations(icases) = length(data_plot.function);
max_case(icases) = max(data_plot.function);
min_case(icases) = min(data_plot.function);
end



fig = randi([0 1000],1);
figure1 = figure(fig);
%figure2 = figure(fig+1);

h = zeros(max(iterations),length(Name)+1);
h2 = zeros(max(iterations),length(Name)+1);


figure(figure1)

set(gca, 'FontName', 'Arial')
set(gca,'FontSize', 18)
set(figure1, 'Position', [0 0 1200 1350])
xlabel('Iterations')
ylabel('Cost Function')

%xlhand = get(gca,'xlabel');
%set(xlhand,'string','X','fontsize',20)
%axis([0 1.05*max(iterations) 0*min(min_case) 1.05*max(max_case)])
axis([0 1.05*max(iterations) 0.95*min(min_case) 1.05*max(max_case)])
axis([0 300 1 5])
%axis([0 1.05*max(iterations) 0.0*min(min_case) 1.05*cost_classical])



classical_valuesx = [0 1.05*max(iterations)];
classical_valuesy = [cost_classical cost_classical];

print(figure1,[Path,Case,'/initial_compliance'],'-dpng')
hold on


for icases = 1:length(Name)
%line_style_plot = line_style{icases};
data_plot.color_plot = data_plot.colors(icases);
data_plot.function = Data{icases}.Cost;

interval = [1:length(data_plot.function)];
h(1:length(interval),icases) = plot_case(icases,Path,Case,interval,data_plot,tipo_selec_micro{icases},figure1,line_style1,classical_valuesx,classical_valuesy);

cost(icases) = data_plot.function(end);
reduction(icases) = abs(cost(icases)-cost(1))/cost(1);

if icases > 1
figure100 = figure(100);
hold on
bar(icases-1,reduction(icases),data_plot.colors(icases))
axis([0 3 0 0.4])
set(gca,'Xtick',1:icases-1,'XTickLabel',data_plot.leyenda2,'FontSize',23,'FontWeight','bold')
set(gca,'YTick',[])
text(icases-1,reduction(icases),[num2str(100*reduction(icases),'%0.1f'),'%'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color',data_plot.colors(icases),'FontSize',30,'FontWeight','bold')
print(figure100,[Path,Case,'/BarReduction',tipo_selec_micro{icases}],'-dpng')
end

end


legend(h(2,1:3),data_plot.leyenda,'FontSize',18,'FontWeight','bold','Location','North')
print(figure1,[Path,Case,'/F_obj'],'-dpng')



 % color string

end

function Data = Obtaindata(Path)
unix(['cd ',Path]);
fid = fopen([Path,'/executed_cases2.txt']);
C = textscan(fid, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
fclose(fid);
k = 1;
for i = 2:2:26
Datos{k} = [C{i}];
k = k + 1;
end
%NewData = reshape(Datos,[],5);
Data.Type = C{16};
Data.Cost = C{8};

end


function h = plot_case(icases,Path,Case,interval,data_plot,name_plot,figure_number,line_style,classical_valuesx,classical_valuesy)

figure(figure_number)
plot(0,data_plot.function(1),'LineWidth',data_plot.linewidth,'Color',data_plot.color_plot,'Marker',data_plot.marker{icases});  
%plot(classical_valuesx,classical_valuesy,'LineWidth',data_plot.linewidth,'Color','g','Marker','*');  
print(figure_number,[Path,Case,'/',name_plot,'_',num2str(1)],'-dpng')

iter_max = length(interval)-1;
iterations_to_print = round(iter_max*((linspace(1,iter_max,7)/iter_max).^2))+1;
iterations_to_print(end+1) = 0;
kprint = 2;

for iplot = 2:length(interval)
intervalo = [iplot-1 iplot];
h(iplot,1) = plot(intervalo-1,data_plot.function(interval(intervalo)));  
line_style_plot = line_style{mod(iplot+1,2)+1};
set(h(iplot,1),'LineWidth',data_plot.linewidth,'Color',data_plot.color_plot,'Marker',data_plot.marker{icases},'LineStyle',line_style_plot)

if iplot-1 == iterations_to_print(kprint)
print(figure_number,[Path,Case,'/',name_plot,'_',num2str(kprint)],'-dpng')
kprint = kprint+1;
end

end

print(figure_number,[Path,Case,'/',name_plot,'_',num2str(kprint)],'-dpng')
end




