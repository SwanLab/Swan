function [ ]=gid_write_headerpost(fid,gtype,ngaus,job)
% write header requered for the post
% 
fprintf(fid,'Gid Post Results File 1.0 \n');
fprintf(fid,'### \n');
fprintf(fid,'# MAT_FEM  V.1.0 \n');
fprintf(fid,'# \n');

fprintf(fid,['GaussPoints "My Gauss" ElemType %s \n'],gtype);
fprintf(fid,['Number Of Gauss Points: %6.0f\n'],ngaus);
if (job==1)
    fprintf(fid,['Natural Coordinates: Internal\n']);
elseif (job==2)
    fprintf(fid,['Natural Coordinates: Given\n']);
    fprintf(fid,'%12.5d %12.5d \n',2/3,1/6);
    fprintf(fid,'%12.5d %12.5d \n',1/6,2/3);
    fprintf(fid,'%12.5d %12.5d \n',1/6,1/6);
end
fprintf(fid,['End gausspoints\n']);

end

