function plotVademecumFull2_2Gid

Path = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Cases/';
Path_2print = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/';
nlevel = 7;
pint = 0; % 0 == no pintar, 1 == pintar
[coordinatesa,phi_snap,theta_snap,numsna] = vademecum_spehere_mesh(nlevel,pint,'ELASTIC');

[conectivities,coordgid] = call_conectivities;
nsnapshots = numsna(end);

for isnap = 1:nsnapshots
    isnap
    fid = fopen([Path,'/',num2str(isnap),'/executed_cases.txt'],'r');     %# Open the file as a binary
    lastline = call_one_before_last_line(fid);
    fclose(fid);  %# Close the file
    C = textscan(lastline, ['%s%s',repmat('%s%f',1,13)]);%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
    k = 1;
    for i = 2:2:28
        Datos{k} = [C{i}];
        k = k + 1;
    end
    type = [Datos{1}];
    Phi = Datos{2};
    Theta = Datos{3};  
    Vol = Datos{4};
    Lagrange0 = Datos{5};
    Iter = Datos{6};
    LagrangeInf = Datos{7};
    Cost0 = Datos{8};
    Cost_Inf = Datos{9};
    h_inf = Datos{10};
    theta_0 = Datos{11};
    theta_inf = Datos{12};
    Gort = Datos{13};
    Mesh = Datos{14};
    
    thet(isnap) = theta_inf;
    Volumen(isnap) = Vol - 0.6;
    ener(isnap) = h_inf;
    lambda_inf(isnap) = LagrangeInf;
    meshes(isnap,1) = Mesh;
    %phi_snap(isnap) = Phi;
    %theta_snap(isnap) = Theta;
    %coordinatesa(isnap,:) = [cos(Phi)*sin(Theta),sin(Phi)*sin(Theta),cos(Theta)];
    iteracio(isnap) = Iter;
    strain(isnap,:) = phitheta2strain(Phi,Theta);
    compact_stress = compact_quadratic_stress(strain(isnap,:));
    
    w1(isnap,1) = compact_stress(1);
    w2(isnap,1) = compact_stress(2);
    w3(isnap,1) = compact_stress(3);
    w4(isnap,1) = compact_stress(4);
    w5(isnap,1) = compact_stress(5);
    w6(isnap,1) = compact_stress(6);
    
     
    fid = fopen([Path,'/',num2str(isnap),'/Ch.txt'],'r');     %# Open the file as a binary
    lastline = call_last_line(fid);
    fclose(fid);  %# Close the file
    C = textscan(lastline, repmat('%s%f',1,11));
    k = 1;
    for i = 2:2:14
        Datos{k} = [C{i}];
        k = k + 1;
    end
    C11(isnap,1) = Datos{2};
    C12(isnap,1) = Datos{3};
    C13(isnap,1) = Datos{4};
    C22(isnap,1) = Datos{5};
    C23(isnap,1) = Datos{6};
    C33(isnap,1) = Datos{7};
    
    determ(isnap,1) = det([C11(isnap,1) C12(isnap,1) C13(isnap,1); C12(isnap,1) C22(isnap,1) C23(isnap,1); C13(isnap,1) C23(isnap,1) C33(isnap,1)]);
 
 
    wC11(isnap,1) = strain(isnap,1)*strain(isnap,1)*C11(isnap,1);
    wC12(isnap,1) = 2*strain(isnap,1)*strain(isnap,2)*C12(isnap,1);
    wC13(isnap,1) = 2*strain(isnap,1)*strain(isnap,3)*C13(isnap,1);
    wC22(isnap,1) = strain(isnap,2)*strain(isnap,2)*C22(isnap,1);
    wC23(isnap,1) = 2*strain(isnap,2)*strain(isnap,3)*C23(isnap,1);
    wC33(isnap,1) = strain(isnap,3)*strain(isnap,3)*C33(isnap,1);
    
    la = coordinatesa(isnap,1);
    lb = coordinatesa(isnap,2);
    lc = coordinatesa(isnap,3);
    lambda1(isnap,1) = real(0.5*((la + lb) - sqrt((la + lb)^2 - 4*(la*lb -lc^2))));
    lambda2(isnap,1) = real(0.5*((la + lb) + sqrt((la + lb)^2 - 4*(la*lb -lc^2))));
    txi(isnap,1) = lambda1(isnap,1)-lambda2(isnap,1);
    txi2(isnap,1) = (lambda1(isnap,1)+lambda2(isnap,1))/2;
    txi3(isnap,1) = lambda1(isnap,1)/lambda2(isnap,1);
    %value = 4*((tan(theta_snap(isnap,1)))*(cos(phi_snap(isnap,1)) + sin(phi_snap(isnap,1))))^2;
    %txi(isnap,1) = (-1 + sqrt(1+value))/(-1 - sqrt(1+value));
    
    fid = fopen([Path_2print,'/vademecum_gid_coordinates.txt'],'a');
    fprintf(fid,'%f %f %f \n',coordinatesa(isnap,1),coordinatesa(isnap,2),coordinatesa(isnap,3));
    fclose(fid);
    
    fid = fopen([Path_2print,'Ch_vademecum.txt'],'a+');
    fprintf(fid,'Isnap %f C11 %f C12 %f C13 %f C22 %f C23 %f C33 %f \n',isnap,C11(isnap,1),C12(isnap,1),C13(isnap,1),C22(isnap,1),C23(isnap,1),C33(isnap,1));
    fclose(fid);
    
    fid = fopen([Path_2print,'wc_vademecum.txt'],'a+');
    fprintf(fid,'Isnap %f w1 %f w2 %f w3 %f w4 %f w5 %f w6 %f \n',isnap,w1(isnap,1),w2(isnap,1),w3(isnap,1),w4(isnap,1),w5(isnap,1),w6(isnap,1));
    fclose(fid);

end
ngaus = 1;   
nnode = 3;
istep = 2;



%[coordinatesa,conevtivities,C11,C12,C13,C22,C23,C33] = 


[phi_full,theta_full] = s1mapping(coordgid(:,2:end)');

%cuadrante1
index_cuad1 = (phi_full <= (pi/4+1e-2) | phi_full >= 2*pi-(pi/4+1e-2)) & theta_full>= 1e-2;



file_name = [Path_2print,'VademecumFull_sphere'];
ToGID_Vademecum(file_name,istep,coordinatesa,conectivities,nnode)
ToGiD_post_Vademecum(file_name,istep,ngaus,ener',determ,C11,C12,C13,C22,C23,C33,wC11,wC12,wC13,wC22,wC23,wC33,w1,w2,w3,w4,w5,w6,txi,lambda1,lambda2,txi2,txi3,lambda_inf,thet,Volumen,meshes,iteracio)

  
end

function [conectivities,coordgid] = call_conectivities
eval('VademecumMeshFull')
conectivities = conectivities(:,2:end);
coordgid = coordinates;
end




function strain = phitheta2strain(phi,theta)

stra_theta = theta;
stra_phi = phi;
stra_x = sin(stra_theta).*cos(stra_phi);
stra_y = sin(stra_theta).*sin(stra_phi);
stra_xy = cos(stra_theta);
strain = [stra_x stra_y stra_xy];

end

