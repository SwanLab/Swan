function [ ] = gid_write_gauss_tensorfield(fid,nameres,time,tfield,idxgp,ngaus,nstre,nelem)
% tfield(ngaus,nstre,nelem)
% i:1..ngaus GID nomenclature
% I:1..ngaus LABFEM nomenclature
% I(i) = idxgp(i)

str = zeros(1,nstre);
% header
if (nstre==3)
s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' Matrix ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
elseif (nstre==4)
    s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' PlainDeformationMatrix ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
else
    fprintf('ERROR in gid_write_gauss_tensorfield. nstre not correctly defined')
end
% components names
fprintf(fid,s,time);
s =['ComponentNames ' '"Sx", ' '"Sy", ' '"Sz", ' '"Sxy", '  '"Syz", ' '"Sxz" ' '\n'];
fprintf(fid,s);

% list of values
fprintf(fid,['Values \n']);
for i=1:nelem
    for istre=1:nstre
        if (ngaus==1)
            str(istre) = tfield(istre,i);
        else
            str(istre) = tfield(idxgp(1),istre,i);
        end
    end
    if (nstre==3)
        fprintf(fid,'%6.0f %12.5d %12.5d %12.5d \n',i,str(:) );
        for ig=2:ngaus
            str = squeeze(tfield(idxgp(ig),:,i));
            fprintf(fid,'%12.5d %12.5d %12.5d \n',str(:) );
        end
    elseif (nstre==4)
        fprintf(fid,'%6.0f %12.5d %12.5d %12.5d %12.5d \n',i,str(:) );
        for ig=2:ngaus
            str = squeeze(tfield(idxgp(ig),:,i));
            fprintf(fid,'%12.5d %12.5d %12.5d %12.5d \n',str(:) );
        end
    end
end
fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

