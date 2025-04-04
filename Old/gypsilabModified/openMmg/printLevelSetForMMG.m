function printLevelSetForMMG(fileName,ls)
fid = fopen(fileName,'w');

%HEADER
fprintf(fid,'%s\n\n','MeshVersionFormatted 2');
fprintf(fid,'%s\n\n',['Dimension ',num2str(dim)]);

% Write nodes
fprintf(fid,'SolAtVertices\n');
fprintf(fid,'[1 1]\n');
fprintf(fid,'\n');

N = size(ls,1);
fprintf(fid,'%d\n',N);
fprintf(fid,'%f %d\n',[ls(:,1),zeros(N,1)]');
fprintf(fid,'\n');

fprintf(fid,'End \n');
end