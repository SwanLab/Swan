function [ ] = gid_write_gauss_vfield(fid,nameres,time,vfield,idxgp)

npnod = size(vfield,1);
s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' Vector ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
fprintf(fid,s,time);

nelem = size(vfield,3);
ngaus = size(vfield,1);

s =['ComponentNames ' '"Sx", ' '"Sy", '  '\n'];
fprintf(fid,['Values \n']);
for i=1:nelem
    fprintf(fid,['%6.0i %13.5d %13.5d \n'],i,squeeze(vfield(idxgp(1),:,i)));
    for ig=2:ngaus
        fprintf(fid,['%13.5d %13.5d \n'],squeeze(vfield(idxgp(ig),:,i)));
    end
end

fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

