function ToGiD_post_Vademecum(file_name,istep,ngaus,Ener,determ,C11,C12,C13,C22,C23,C33,wC11,wC12,wC13,wC22,wC23,wC33,w1,w2,w3,w4,w5,w6,txi,lambda1,lambda2,txi2,txi3,lambda_inf,theta_res,Volumen,mesh,iteracio)

gtype = 'Triangle'; %gid type

% Escribe el fichero de resultados
res_file = strcat(file_name,'_',num2str(istep),'.flavia.res');
fid = fopen(res_file,'w');


job=2;
gid_write_headerpost(fid,gtype,ngaus,job)

timestep = 1;

nameres = 'Energia';
gid_write_nodal_sclfield(fid,nameres,timestep,Ener);

nameres = 'Determinante Ch';
gid_write_nodal_sclfield(fid,nameres,timestep,determ);

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

nameres = 'wC11';
gid_write_nodal_sclfield(fid,nameres,timestep,wC11);

nameres = 'wC12';
gid_write_nodal_sclfield(fid,nameres,timestep,wC12);

nameres = 'wC13';
gid_write_nodal_sclfield(fid,nameres,timestep,wC13);

nameres = 'wC22';
gid_write_nodal_sclfield(fid,nameres,timestep,wC22);

nameres = 'wC23';
gid_write_nodal_sclfield(fid,nameres,timestep,wC23);

nameres = 'wC33';
gid_write_nodal_sclfield(fid,nameres,timestep,wC33);



nameres = 'w1';
gid_write_nodal_sclfield(fid,nameres,timestep,w1);

nameres = 'w2';
gid_write_nodal_sclfield(fid,nameres,timestep,w2);

nameres = 'w3';
gid_write_nodal_sclfield(fid,nameres,timestep,w3);

nameres = 'w4';
gid_write_nodal_sclfield(fid,nameres,timestep,w4);

nameres = 'w5';
gid_write_nodal_sclfield(fid,nameres,timestep,w5);

nameres = 'w6';
gid_write_nodal_sclfield(fid,nameres,timestep,w6);




nameres = 'txi';
gid_write_nodal_sclfield(fid,nameres,timestep,txi);

nameres = 'lambda1';
gid_write_nodal_sclfield(fid,nameres,timestep,lambda1);

nameres = 'lambda2';
gid_write_nodal_sclfield(fid,nameres,timestep,lambda2);

nameres = 'txi2';
gid_write_nodal_sclfield(fid,nameres,timestep,txi2);

nameres = 'txi3';
gid_write_nodal_sclfield(fid,nameres,timestep,txi3);

nameres = 'lambda_inf';
gid_write_nodal_sclfield(fid,nameres,timestep,lambda_inf');

nameres = 'theta_res';
gid_write_nodal_sclfield(fid,nameres,timestep,theta_res');

nameres = 'Volumen';
gid_write_nodal_sclfield(fid,nameres,timestep,Volumen');

nameres = 'mesh';
gid_write_nodal_sclfield(fid,nameres,timestep,mesh);

nameres = 'Number of Iterations';
gid_write_nodal_sclfield(fid,nameres,timestep,iteracio);


status = fclose(fid);

end
