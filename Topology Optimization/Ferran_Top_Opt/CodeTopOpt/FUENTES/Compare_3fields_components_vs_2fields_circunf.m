function Compare_3fields_components_vs_2fields_circunf


Case_3Fields = 'Airfoil_ArribaIzayDerEmp_abajoder_hor_fij_5grupos_micro_hor_punto_gauss';
Case_2Fields = 'Airfoil_ArribaIzayDerEmp_abajoder_hor_fij_5grupos_micro_hor';

%Case_3Fields = 'Airfoil_ArribaIzayDerEmp_abajoder_hor_fij_5grupos_punto_gauss';
%Case_2Fields = 'Airfoil_ArribaIzayDerEmp_abajoder_hor_fij_5grupos';
   
%Case_3Fields = 'CantilibearBeam3Fields_and_components';%'Airfoil3Fields_and_components/';%'Cantilibear3Fields_and_components/';
%Case_2Fields = 'CantilibearBeam2Fields_Circunf';%'Airfoil2Fields_Circunf/';%'Cantilibear2Fields_Circunf/';

%Case_3Fields = 'Airfoil3Fields_and_components/';%'Cantilibear3Fields_and_components/';
%Case_2Fields = 'Airfoil2Fields_Circunf/';%'Cantilibear2Fields_Circunf/';
Path = '/home/aferrer/Documents/Doctorat/Tesi/FEM_DT/Results/';

%Case_3Fields = 'Lambda0_pen_1_factor_0_99_V_0_6_3Fields';%'Airfoil3Fields_and_components/';%'Cantilibear3Fields_and_components/';
%Case_2Fields = 'Lambda0_pen_1_factor_0_99_V_0_6_2Fields';%'Airfoil2Fields_Circunf/';%'Cantilibear2Fields_Circunf/';
%Path = '/home/aferrer/Desktop/PreubasAirfoil/';


Case_1 = Case_3Fields;
Case_2 = Case_2Fields;

%Case_1 = 'ArribaIzayDerEmp_abajoder_hor_fij_5grupos_micro_hor';
%Case_2 = 'ArribaIzayDerEmp_abajoder_hor_fij_5grupos_micro_hor2';
%Case_3 = 'ArribaIzayDerEmp_abajoder_hor_fij_5grupos_micro_hor3';


%Path = '/home/aferrer/Desktop/PreubasAirfoil/';


Name =  {Case_1,Case_2};
%Name =  {Case_1,Case_2,Case_3};
leyenda = {'3 Fileds and components','2 Fields with circle'};
%leyenda = {'penalty 1, 0.9 ','penalty 10, 0.99','penalty 1, 0.99'};
colors = 'brgcmyk';
ending = [83 58];

fig = randi([0 1000],1);
figure(fig)

for icases = 1:length(Name)
Cost = Obtaindata([Path,Name{icases}]);    

h(icases) = plot(Cost(1:ending(icases)),'+-');

set(h(icases),'LineWidth',2,'Color',colors(icases))

hold on
end
legend(leyenda,'FontSize',8,'FontWeight','bold')
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

end




