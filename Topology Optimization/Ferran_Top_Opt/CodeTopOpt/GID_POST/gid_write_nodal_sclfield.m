function [ ] = gid_write_nodal_sclfield(fid,nameres,time,sclfield)

npnod = size(sclfield,1);
s =['Result' ' "' nameres '" ' '"time"' '%12.5d' ' Scalar ' ' OnNodes ' '\n'];
fprintf(fid,s,time);
fprintf(fid,['Values \n']);
% for i = 1 : npnod
%    fprintf(fid,['%6.0i %12.5d \n'],i,sclfield(i));
% end
fprintf(fid,['%6.0i %12.5d \n'],[(1:npnod);sclfield']); % vectorized NEW! (Ferran)

fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end


