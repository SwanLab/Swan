function plotVademecum2Gid
Path = '/media/Alex/Vademecum';
nlevel = 7;
pint = 0; % 0 == no pintar, 1 == pintar
[coordinatesa,phi_snap,theta_snap,numsna] = vademecum_spehere_mesh(nlevel,pint);

[conectivities,coordgid] = call_conectivities;
nsnapshots = numsna(end);

for isnap = 1:nsnapshots
    isnap
    fid = fopen([Path,'/',num2str(isnap),'/executed_cases.txt'],'r');     %# Open the file as a binary
    lastline = call_last_line(fid);
    fclose(fid);  %# Close the file
    C = textscan(lastline, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
    k = 1;
    for i = 2:2:28
        Datos{k} = [C{i}];
        k = k + 1;
    end
    tipo = [Datos{1}];
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
    Energia = Datos{13};
    Mesh = Datos{14};
    
    thet(isnap) = theta_inf;
    Volumen(isnap) = Vol - 0.6;
    ener(isnap) = Energia;
    lambda_inf(isnap) = LagrangeInf;
    meshes(isnap,1) = Mesh;
    %phi_snap(isnap) = Phi;
    %theta_snap(isnap) = Theta;
    %coordinatesa(isnap,:) = [cos(Phi)*sin(Theta),sin(Phi)*sin(Theta),cos(Theta)];
    iteracio(isnap) = isnap;
    
    
    fid = fopen([Path,'/',num2str(isnap),'/Ch.txt'],'r');     %# Open the file as a binary
    lastline = call_last_line(fid);
    fclose(fid);  %# Close the file
    C = textscan(lastline, '%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
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
    
    fid = fopen([Path,'/vademecum_gid_coordinates.txt'],'a');
    fprintf(fid,'%f %f %f \n',coordinatesa(isnap,1),coordinatesa(isnap,2),coordinatesa(isnap,3));
    fclose(fid);

end
ngaus = 1;   
nnode = 3;
istep = 2;

%[coordinatesa,conevtivities,C11,C12,C13,C22,C23,C33] = 

file_name = 'Vademecum';
ToGID_Vademecum(file_name,istep,coordinatesa,conectivities,nnode)
ToGiD_post_Vademecum(file_name,istep,ngaus,ener',C11,C12,C13,C22,C23,C33,txi,lambda1,lambda2,txi2,txi3,lambda_inf,thet,Volumen,meshes)



end

function [conectivities,coordgid] = call_conectivities
eval('VademecumMesh')
conectivities = conectivities(:,2:end);
coordgid = coordinates;
end


function lastLine = call_last_line(fid)

lastLine = '';                   %# Initialize to empty
offset = 1;                      %# Offset from the end of file
fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
newChar = fread(fid,1,'*char');  %# Read one character
while (~strcmp(newChar,char(10))) || (offset == 1)
    lastLine = [newChar lastLine];   %# Add the character to a string
    offset = offset+1;
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    newChar = fread(fid,1,'*char');  %# Read one character
end

end