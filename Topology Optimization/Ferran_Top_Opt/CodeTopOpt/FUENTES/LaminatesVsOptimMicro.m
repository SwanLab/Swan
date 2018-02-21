function LaminatesVsOptimMicro

physical_type = 'THERMAL';%'ELASTIC';%'THERMAL''ELASTIC';%'THERMAL'
pint = 0; % 0 == no pintar, 1 == pintar
nlevel = 7;
[def,phi_all,theta_all,numsnapshots] = vademecum_spehere_mesh(nlevel,pint,physical_type);



printi = 0;

[~,index_sort] = sort(phi_all);


for iphi = 1:length(phi_all)
phi_date = phi_all(iphi);    
[Ch,fobj] = obtainCh_Fobj('optimo',phi_date,iphi);
Ch_opt(:,:,iphi) = Ch;
fobj_opt(iphi,1) = fobj;

[Ch,fobj] = obtainCh_Fobj('LaminadoNoPeriodico',phi_date,iphi);
Ch_LamNo(:,:,iphi) = Ch;
fobj_LamNo(iphi,1) = fobj;

[Ch,fobj] = obtainCh_FobjRotated(phi_date);
Ch_Lam(:,:,iphi) = Ch;
fobj_Lam(iphi,1) = fobj;
%[matCh,h0] = evaluate_laminates(phi_date,printi,physical_type,iphi);
end



plotfunction(phi_all(index_sort),fobj_opt(index_sort),fobj_LamNo(index_sort),fobj_Lam(index_sort),'Cost Function','Fobj_laminatesVsOptim')
plotfunction(phi_all(index_sort),squeeze(Ch_opt(1,1,index_sort)), squeeze(Ch_LamNo(1,1,index_sort)), squeeze(Ch_Lam(1,1,index_sort)),'K11','K11_laminatesVsOptim')
plotfunction(phi_all(index_sort),squeeze(Ch_opt(1,2,index_sort)), squeeze(Ch_LamNo(1,2,index_sort)), squeeze(Ch_Lam(1,2,index_sort)),'K12','K12_laminatesVsOptim')
plotfunction(phi_all(index_sort),squeeze(Ch_opt(2,2,index_sort)), squeeze(Ch_LamNo(2,2,index_sort)), squeeze(Ch_Lam(2,2,index_sort)),'K22','K22_laminatesVsOptim')

end

function plotfunction(x,y1,y2,y3,title,nameFile)
leyenda = {'Optima DT','Laminado No periodico','Laminado'};
h = figure(1);
plot(x,[y1 y2 y3],'-+')
figure(h)

set(gca, 'FontName', 'Arial')
set(gca,'FontSize', 15)
set(h, 'Position', [0 0 1200 1350])
xlabel('Angle')
ylabel(title)

%set(h(icases),'LineWidth',2,'Color',colors(icases),'LineStyle','-+')
legend(leyenda,'FontSize',8,'FontWeight','bold','Location','northwest')
print(h,['/home/aferrer/Documents/Doctorat/Tesi/Rastro/figures/',nameFile],'-dpng')



end




function changeFirstLineCh(Path)
unix(['cd ',Path]);
fid = fopen([Path,'/Ch.txt']);
tline = fgetl(fid);
fclose(fid);
fid = fopen([Path,'/Ch.txt'],'r+');
linewrite = [tline(1:8), num2str(0), tline(9:end)];
fprintf(fid,repmat('%s%f',1,8),linewrite);
fprintf(fid,'\n');
fclose(fid);



end

function [Ch,fobj] = obtainCh_FobjRotated(phi_date)
alpha = unitary_macro_flux_vector(phi_date);
R = [cos(phi_date) -sin(phi_date);sin(phi_date) cos(phi_date)]';
file_path = ['/home/aferrer/Documents/Doctorat/Tesi/ElasticThermal_Regularized','/VademecumThermal_','optimo','/',num2str(2),'/'];
Ch_hor = ObtaindataCh(file_path);
Ch = R'*Ch_hor*R;
fobj = objectivefunction(alpha,Ch);

end

function [Ch,fobj] = obtainCh_Fobj(case_study,phi_date,iphi)
file_path = ['/home/aferrer/Documents/Doctorat/Tesi/ElasticThermal_Regularized','/VademecumThermal_',case_study,'/',num2str(iphi),'/'];
alpha = unitary_macro_flux_vector(phi_date);
Ch = ObtaindataCh(file_path);
fobj = objectivefunction(alpha,Ch);

end

function matCh = ObtaindataCh(Path)
unix(['cd ',Path]);
fid = fopen([Path,'/Ch.txt']);
C = textscan(fid,repmat('%s%f',1,8));
fclose(fid);
k = 1;
for i = 4:2:15
Datos{k} = [C{i}];
k = k + 1;
end
%NewData = reshape(Datos,[],5);
matCh =  [Datos{1}(end) Datos{2}(end); Datos{2}(end) Datos{3}(end)];
%Cost = Datos{8};

end

function alpha = unitary_macro_flux_vector(phi_date)
alpha = [cos(phi_date) sin(phi_date)]';
end
function fobj = objectivefunction(alpha,Ch)
fobj = alpha'*(Ch\alpha);
end
