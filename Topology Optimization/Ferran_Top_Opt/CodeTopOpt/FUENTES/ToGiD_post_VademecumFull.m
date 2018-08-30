function ToGiD_post_VademecumFull(file_name,istep,ngaus,Ener,C11,C12,C13,C22,C23,C33)

gtype = 'Triangle'; %gid type

% Escribe el fichero de resultados
res_file = strcat(file_name,'_',num2str(istep),'.flavia.res');
fid = fopen(res_file,'w');


job=2;
gid_write_headerpost(fid,gtype,ngaus,job)

timestep = 1;

nameres = 'Energia';
gid_write_nodal_sclfield(fid,nameres,timestep,Ener);

nameres = 'C11';
gid_write_nodal_sclfield(fid,nameres,timestep,C11);

nameres = 'C12';
gid_write_nodal_sclfield(fid,nameres,timestep,C12);

nameres = 'C13';
gid_write_nodal_sclfield(fid,nameres,timestep,C13);

nameres = 'C22';
gid_write_nodal_sclfield(fid,nameres,timestep,C22);

nameres = 'C23';
gid_write_nodal_sclfield(fid,nameres,timestep,C23);

nameres = 'C33';
gid_write_nodal_sclfield(fid,nameres,timestep,C33);

status = fclose(fid);

end
