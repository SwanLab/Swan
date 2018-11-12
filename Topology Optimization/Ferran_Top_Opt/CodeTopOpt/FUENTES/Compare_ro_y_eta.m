function Compare_ro_y_eta


Path = '/home/aferrer/Desktop/PreubasAirfoil/';


colors = 'brgcmyk';
Marker = {'-','--'};


eta = [0.99 0.999];
ro = [ 1 5 10];
leyenda = {'ro = 1','ro = 5','ro = 10'};


for ieta = 1:length(eta)


fig = randi([0 1000],1);
figure_number = figure(fig);


for icases = 1:length(ro)
    
total_path = [Path,'Cantilever_ro',strrep(num2str(ro(icases)),'.','_'),'_eta',strrep(num2str(eta(ieta)),'.','_'),'/'];    
Cost = Obtaindata(total_path);    

%h(icases) = plot(Cost(1:ending(icases)),'+-');
h(icases) = plot(Cost);

set(h(icases),'LineWidth',2,'Color',colors(icases),'LineStyle',Marker{ieta})
legend(leyenda,'FontSize',8,'FontWeight','bold')

set(gca, 'FontName', 'Arial')
set(gca,'FontSize', 18)
xlabel('Iterations')
ylabel('Cost Function')


hold on
end
name_plot = ['F_obj_eta',strrep(num2str(eta(ieta)),'.','_'),'diferent_ro_reduced'];
print(figure_number,[Path,'/',name_plot],'-dpng')


end

%legend([';''])


 % color string

end

function Cost = Obtaindata(Path)
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
Cost = Datos{4};
%Cost = Datos{8};

end




