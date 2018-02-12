function [Ch,phi,theta,phi_matrix,theta_matrix,number_matrix,theta_step] = call_Ch(phyisical_type)
%fid = fopen('/home/aferrer/Documents/Doctorat/Tesi/MicroEscalaDT/Ch_vademecum.txt','r'); %# Open the file as a binary
fid = fopen('/home/aferrer/Documents/Doctorat/Tesi/MicroStructures_Vademecum/Vademecum_Penalty1_Thetamin1/Ch_vademecumQuarter.txt','r'); %# Open the file as a binary

C = textscan(fid,repmat('%s%f',1,7));
fclose(fid);  %# Close the file
Ch = C{4:2:14};
C11 = C{4};
C12 = C{6};
C13 = C{8};
C22 = C{10};
C23 = C{12};
C33 = C{14};
Ch = [C11 C12 C13 C22 C23 C33];

pint = 0; nlevel = 7; % 0 == no pintar, 1 == pintar
[~,phi,theta,numsnapshots,phi_matrix,theta_matrix,number_matrix,theta_step] = vademecum_spehere_mesh(nlevel,pint,phyisical_type);
end