function [ ] = gid_write_tensorfield(fid,nameres,time,tfield)

nelem = size(tfield,1);
nstre = size(tfield,2);
s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' Matrix ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
fprintf(fid,s,time);
s =['ComponentNames ' '"Sx", ' '"Sy", ' '"Sxy", ' '"Sz", '  '"Syz", ' '"Sxz" ' '\n'];
fprintf(fid,s);

fprintf(fid,['Values \n']);
for i = 1 : nelem
   str(1) = tfield(1,i);
   str(2) = tfield(2,i);
   str(3) = tfield(3,i);
   str(4) = tfield(4,i);
   fprintf(fid,'%6.0f %12.5d %12.5d  %12.5d %12.5d  0.0  0.0 \n',i,str(:) );
end
fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

