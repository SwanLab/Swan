function [conectivities,coordgid] = create_whole_mesh_shpere
[conectivities,coordgid] = call_conectivities_1_over_8_sphere;
Path = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Cases/';
Path_2print = '/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/';
nSnap = size(coordgid,1);


flag_read_Ch = 1;
flag_read_all_variables = 1;
flag_printxyzpoint = 0;

if flag_read_all_variables == 1
    [ener,phi,theta,strain] = readAllData(Path,nSnap,coordgid);
else
    load_var = load('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/AllVariablesVademecum');
    ener = load_var.ener;
    phi = load_var.Phi;
    theta = load_var.Theta;
    strain = load_var.strain;
end




if flag_read_Ch == 1
    Ch = readCh(Path,nSnap);
else
    loadCh = load('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Ch');
    Ch = loadCh.Ch;
end

%[phi,theta] = s1mapping(coordgid(:,2:end)');
%theta = mod(theta+pi,pi);
%index_border = (phi >= pi/4 - 1e-3 & phi <= pi/4 + 1e-3) | (phi <= (7*pi/4 + 1e-3) & phi >= 7*pi/4 - 1e-3)  | (theta <= 1e-3 | theta >= pi/2 -1e-3);
%index_interior = 1-index_border;   

[cuadrante,coord_full,Ch_full,ener_full] = obtainFullVademecumVariables(coordgid,Ch,ener,phi,theta,strain);

if flag_printxyzpoint == 1
printxyzpoints(Path_2print,coord_full)
end


ngaus = 1;   
nnode = 3;
istep = 2;
[conectivities_full,coordgid_full] = call_conectivities_full;

C11 = squeeze(Ch_full(1,1,:));
C12 = squeeze(Ch_full(1,2,:));
C13 = squeeze(Ch_full(1,3,:));
C22 = squeeze(Ch_full(2,2,:));
C23 = squeeze(Ch_full(2,3,:));
C33 = squeeze(Ch_full(3,3,:));



file_name = [Path_2print,'VademecumFull_sphere'];
ToGID_Vademecum(file_name,istep,coordgid_full(:,2:end),conectivities_full,nnode)
ToGiD_post_VademecumFull(file_name,istep,ngaus,ener_full',C11,C12,C13,C22,C23,C33)

end

function [conectivities_full,coordgid_full] = call_conectivities_full
eval('VademecumMeshFull')
conectivities_full = conectivities(:,2:end);
%conectivities_full(:,1) = 1:size(conectivities_full,1);
coordgid_full = coordinates;

end

function printxyzpoints(Path_2print,coord_full)
fid = fopen([Path_2print,'/vademecumFull_gid_coordinates.txt'],'w+');
for isnap = 1: size(coord_full,1)
    fprintf(fid,'%f %f %f \n',coord_full(isnap,1),coord_full(isnap,2),coord_full(isnap,3));
end
 fclose(fid);
end

function Ch = readCh(Path,nSnap)
Ch = zeros(3,3,nSnap);
for isnap = 1:nSnap

fid = fopen([Path,'/',num2str(isnap),'/Ch.txt'],'r');     %# Open the file as a binary
    lastline = call_last_line(fid);
    fclose(fid);  %# Close the file
    C = textscan(lastline, repmat('%s%f',1,11));
    k = 1;
    for i = 2:2:14
        Datos{k} = [C{i}];
        k = k + 1;
    end
    C11 = Datos{2};
    C12 = Datos{3};
    C13 = Datos{4};
    C22 = Datos{5};
    C23 = Datos{6};
    C33 = Datos{7};
    
    Ch(:,:,isnap) = [ C11 C12 C13; C12 C22 C23; C13 C23 C33];
    
    fid = fopen([Path,'Ch_vademecum.txt'],'a+');
    fprintf(fid,'Isnap %f C11 %f C12 %f C13 %f C22 %f C23 %f C33 %f \n',isnap,C11,C12,C13,C22,C23,C33);
    fclose(fid);
    
end

save('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Ch');

end



function [ener,Phi,Theta,strain] = readAllData(Path,nSnap,coordinatesa)
for isnap = 1:nSnap
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
    Phi(isnap) = Datos{2};
    Theta(isnap) = Datos{3};  
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
    strain(isnap,:) = phitheta2strain(Phi(isnap),Theta(isnap));
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
    
end

save('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/AllVariablesVademecum');

end



