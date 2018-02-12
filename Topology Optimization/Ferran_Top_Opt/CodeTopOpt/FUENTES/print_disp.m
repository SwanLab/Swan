function [  ] = print_disp(fid,timestep,nunkn,npnod,nbdata,vdisp)
% Imprime los desplazamientos de los tres experimentos

nndof = npnod*nunkn;
for idata=1:nbdata
    d_u = vdisp(idata,:);
    tdisp(1:npnod,1)=d_u(1:nndof);
    nameres = ['DISP EXP ' num2str(idata)];
    gid_write_nodal_sclfield(fid,nameres,timestep,tdisp);
end

