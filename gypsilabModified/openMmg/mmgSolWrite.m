function mmgSolWrite(filename,sol,is2d)
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique
% Algiane Froehly, CARDAMOME, INRIA-SOFT
% LGPL Lesser General Public License v3.0.
% Remeshing using Mmg tools : https://www.mmgtools.org

% Open .sol file
fid = fopen(filename,'w');
if(fid==-1)
    error('mmgSolWrite.m : unavailable case');
end

% Header
fprintf(fid,'%s\n','# SOL generated using Gyspilab for Matlab #');
fprintf(fid,'%s\n\n','MeshVersionFormatted 2');
if is2d
    fprintf(fid,'%s\n\n','Dimension 2');
else
    fprintf(fid,'%s\n\n','Dimension 3');
end

% Solutions at vertices only
fprintf(fid,'SolAtVertices\n');
fprintf(fid,'%d\n',length(sol));

% Scalar : 1
if (size(sol,2)==1)
    fprintf(fid,'1 1\n');
    fprintf(fid,'%12.8f\n',sol');
    
% Vectorial : 2 
elseif (size(sol,2)==3)
    fprintf(fid,'1 2\n');
    fprintf(fid,'%12.8f %12.8f %12.8f\n',sol');

% Tensor : 3 (symetric, lower part)
elseif (size(sol,2)==6)
    fprintf(fid,'1 3\n');
    fprintf(fid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',sol');

% Other    
else
    error('mmgSolWrite.m : unavailable case');
end
fprintf(fid,'\n');

% Termination
fprintf(fid,'End \n');

% Close mesh file
fclose(fid);
end
