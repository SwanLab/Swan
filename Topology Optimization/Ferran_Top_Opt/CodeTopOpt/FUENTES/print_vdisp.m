function [  ] = print_vdisp(fid,timestep,ndime,npnod,nbdata,vdisp)
% Imprime los desplazamientos de los tres experimentos

nndof = npnod*ndime;
for idata=1:nbdata
    d_u = vdisp(idata,:);
    tdisp(1:npnod,1)=d_u(1:2:nndof);
    tdisp(1:npnod,2)=d_u(2:2:nndof);
    nameres = ['DISP EXP ' num2str(idata)];
    gid_write_vfield(fid,nameres,timestep,tdisp);
end

