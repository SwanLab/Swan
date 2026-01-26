function [StwoST,celasST]= SmallStrainLargeRotations(EgreenlST,MATPRO)

StwoST  = zeros(size(EgreenlST)) ;
nstrain = size(MATPRO.celasglo,2) ;
INDEX = cell(1,nstrain) ;
for istrain = 1:nstrain
    INDEX{istrain} = istrain:nstrain:length(StwoST) ;
end
for  i = 1:nstrain
    for j= 1:nstrain
        StwoST(INDEX{i}) = StwoST(INDEX{i}) + MATPRO.celasglo(INDEX{i},j).*EgreenlST(INDEX{j}) ;
    end
end
celasST  = MATPRO.celasglo ;