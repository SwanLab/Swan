function replot_values_Ch

Path_m = ['/media/Alex/Vademecum/'];
Path_p = ['/home/aferrerferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/'];
file_name = 'RVE04N3';
nsnap = 2145;


for isnap = 2043:nsnap
    isnap
Path_mesh = [Path_m,num2str(isnap),'/'];
Path_phi = [Path_p,num2str(isnap)];
nmesh = get_number_mesh(Path_mesh);
phifunct = get_phi(Path_phi);
matCh = get_Ch(file_name,nmesh,phifunct);
printCh(matCh,isnap,Path_mesh)
end













end

function printCh(matCh,isnap,Path)
fid = fopen([Path,'/Ch.txt'],'a+');     %# Open the file as a binary
fprintf(fid,'SNAP: %3.0d  C11 %3.6f C12 %3.6f C13 %3.6f C22 %3.6f C23 %3.6f C33 %3.6f  \n',...
    isnap,matCh(1,1),matCh(1,2),matCh(1,3),matCh(2,2),matCh(2,3),matCh(3,3));
fclose(fid);

fid = fopen('/home/aferrerferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Ch_vademecum','a+');
fprintf(fid,'SNAP: %3.0d  C11 %3.8f C12 %3.8f C13 %3.8f C22 %3.8f C23 %3.8f C33 %3.8f \n',...
    isnap,matCh(1,1),matCh(1,2),matCh(1,3),matCh(2,2),matCh(2,3),matCh(3,3));
fclose(fid);

end

function [fext,element,fixnodes,problembsc,coordinates,phifunct_n] = call_read_data(file_name,lambda,Vfrac,penalty,imesh)

[fext,element,fixnodes,problembsc,coordinates,phifunct_n] = read_data_problem(file_name,imesh);
              
element.material.opt_L = lambda; % bulk 20;
element.material.Vfrac = Vfrac;
element.material.penalty = penalty;
end


function matCh = get_Ch(file_name,nmesh,phifunct_n)
[fext,element,fixnodes,problembsc,coordinates,~] = call_read_data(file_name,0,0.6,4,nmesh);
% Initialize basics dimensions
[dim.npnod,dim.nndof,dim.ndime] = data_nod(coordinates, element.type,problembsc.problemtype);
[dim.nelem,dim.nnode,dim.neleq] = data_elem(element.conectivities, element.type);

% basic variables
% Initial value of 
fext.h_C_0 = 1;
[phigp_n] = interpol(phifunct_n,element,dim,problembsc);
fext.macrostra = [1 0 0];
[~,matCh] = module_M(phifunct_n,phigp_n,element,fixnodes,problembsc,coordinates,fext,dim);
end



function phi = get_phi(Path)
fid = fopen([Path,'.res'], 'rt');
% read the entire file, if not too big
s = textscan(fid, '%s', 'delimiter', '\n');
% search for your Region:
idx1 = find(strcmp(s{1}, 'Values '), 1, 'first');
% now search for the next Region:
idx2 = find(strcmp(s{1}, 'End Values '), 1, 'first');

node = [idx1+1:idx2-1];
for inode = 1:length(node)
S = str2num(s{1}{node(inode)});
phi(inode,1) = S(2);
end
%S2 = s{1}{idx2}
fclose(fid);


end


function nmesh = get_number_mesh(Path)

fid = fopen([Path,'/executed_cases.txt'],'r');     %# Open the file as a binary
lastline = call_last_line(fid);
fclose(fid);  %# Close the file
C = textscan(lastline, '%s%s%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f');
% k = 1;
% for i = 28:2:28
%     Datos{k} = [C{i}];
%     k = k + 1;
% end
% tipo = [Datos{1}];
% Phi = Datos{2};
% Theta = Datos{3};
% Vol = Datos{4};
% Lagrange0 = Datos{5};
% Iter = Datos{6};
% LagrangeInf = Datos{7};
% Cost0 = Datos{8};
% Cost_Inf = Datos{9};
% h_inf = Datos{10};
% theta_0 = Datos{11};
% theta_inf = Datos{12};
% Energia = Datos{13};
Mesh = C{28};
nmesh = Mesh;



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