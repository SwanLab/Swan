function [MASTER SLAVES] =FindLocationMasterSlave(MASTER,SLAVES,COOR)
if nargin == 0
    load('tmp.mat')
end

for i=1:length(MASTER)
    Nmast = MASTER{i} ;
    coorMAST = COOR(Nmast,:) ;  % Master coordinates   
    for j=1:size(SLAVES,2)
        Nslv = SLAVES{i,j} ;
        if ~isempty(Nslv) && length(Nmast)>1
            coorSLV = COOR(Nslv,:) ;  % Slave coordinates
            for inode = 1:length(Nslv)
                
            end
        end
    end
end