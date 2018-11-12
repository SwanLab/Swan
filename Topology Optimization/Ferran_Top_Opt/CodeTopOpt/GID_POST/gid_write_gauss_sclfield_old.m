function [ ] = gid_write_gauss_sclfield(fid,nameres,time,sfield,idxgp)

nelem = size(sfield,2);
ngaus = size(sfield,1);

% header
s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' Scalar ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
fprintf(fid,s,time);
% list of values
fprintf(fid,['Values \n']);
for i=1:nelem
    fprintf(fid,'%6.0f %12.5d \n',i,sfield(idxgp(1),i));
    for ig=2:ngaus
        fprintf(fid,'%12.5d \n',sfield(idxgp(ig),i));
    end
end
fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

