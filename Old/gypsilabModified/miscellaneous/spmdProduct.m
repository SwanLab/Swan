function MV = spmdProduct(Ilab,M,V)
% Copyright (c) 2019, Matthieu Aussal, Ecole Polytechnique, CMAP      
% GNU General Public License v3.0. 
% Matrix vector product using SPMD repartition 

% Parallel Matrix-Vector product
spmd
    if iscell(M)
        MV = 0;
        for j = 1:numlabs
            MV = MV + M{j}*V(Ilab{j},:);
        end
    else
        MV = M*V;
    end
end

% Data recuperation
MV = cell2mat(MV(:));

% Reordering
MV(cell2mat(Ilab),:) = MV;
end