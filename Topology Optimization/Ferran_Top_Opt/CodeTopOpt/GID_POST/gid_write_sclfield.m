function [ ] = gid_write_sclfield(fid,nameres,time,sclfield)

nelem = size(sclfield,1);
s =['Result' ' "' nameres '" ' '"time"' '%12.5d' ' Scalar ' ' OnGaussPoints ' ...
    '"' 'My Gauss' '"' '\n'];
fprintf(fid,s,time);
fprintf(fid,['Values \n']);
for i = 1 : nelem
   fprintf(fid,['%6.0i %12.5d \n'],i,sclfield(i));
end
fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end


