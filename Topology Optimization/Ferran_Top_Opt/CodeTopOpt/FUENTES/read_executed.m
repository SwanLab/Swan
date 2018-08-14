function read_executed

fini = 141.5203;


Path{1} = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/Macro_Ch_vademecum_random_2/';
Path{2} = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/Cantilibear2Fields_Circunf/';
Path{3} = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/Cantilibear3Fields_ini_circunf/';

for i = 1:3
file = [Path{i},'executed_cases.txt'];
fid = fopen(file,'r');
formato = repmat('%s%f',1,10);
i
C = textscan(fid, formato);
fclose(fid);
k = 1;

h{i} = C{8};
Vol{i} = C{14};
end
figure(3)
% plot([fini;h{1}])
% hold on
plot([fini;h{2}],'*-r','LineWidth',3,'MarkerSize',10)
hold on
plot([fini;h{3}],'*-b','LineWidth',3,'MarkerSize',10)
set(gca,'fontsize',35)

figure(2)
% plot([fini;h{1}])
% hold on
plot([1;Vol{2}],'*-r','LineWidth',3,'MarkerSize',10)
hold on
plot([1;Vol{3}],'*-b','LineWidth',3,'MarkerSize',10)
set(gca,'fontsize',35)

end

