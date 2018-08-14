function [ ] = gid_write_gauss_vfield(fid,nameres,time,vfield,idxgp,ndime)

% npnod = size(vfield,1);
s =['Result' ' "' nameres '" ' '"time" ' '%12.5d' ' Vector ' ' OnGaussPoints ' '"' 'My Gauss' '"' '\n'];
fprintf(fid,s,time);

nelem = size(vfield,3);
ngaus = size(vfield,1);

% s =['ComponentNames ' '"Sx", ' '"Sy", '  '\n'];
% fprintf(fid,s);
fprintf(fid,['Values \n']);
% for i=1:nelem
%     fprintf(fid,['%6.0i %13.5d %13.5d \n'],i,squeeze(vfield(idxgp(1),:,i)));
%     for ig=2:ngaus
%         fprintf(fid,['%13.5d %13.5d \n'],squeeze(vfield(idxgp(ig),:,i)));
%     end
% end


%He desactivado esta PARTE
% for i=1:nelem
%     fprintf(fid,['%6.0i '],i);
%     fprintf(fid,['%13.5d %13.5d \n'],squeeze(vfield(idxgp(1:ngaus),:,i))); %fprintf('%13.5d %13.5d \n',permute(a,[2 1 3]))
% %     for ig=2:ngaus
% %         fprintf(fid,['%13.5d %13.5d \n'],squeeze(vfield(idxgp(ig),:,i)));
% %     end
% end
%He desactivado esta PARTE

vfield = vfield(idxgp,:,:);
vfield = reshape(permute(vfield,[2 1 3]),ngaus*ndime,[]);
format = ['%d',repmat([repmat(' %f',1,ndime),'\n',' '],1,ngaus)];
fprintf(fid,format,[1:nelem;vfield]);

fprintf(fid,['End Values \n']);
fprintf(fid,'# \n');

end

