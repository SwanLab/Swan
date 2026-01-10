function [NODES_n,GAUSSV_n] = UpdateHistoricalVariables(NODESV_np1,GAUSSV_np1)

fff  =fieldnames(NODESV_np1);
for  ifield = 1:length(fff)
    FFF = fff{ifield} ;
    NODES_n.(FFF) = NODESV_np1.(FFF) ;
end
fff  =fieldnames(GAUSSV_np1) ; 
for  ifield = 1:length(fff)
    FFF = fff{ifield} ;
    GAUSSV_n.(FFF) = GAUSSV_np1.(FFF) ;
end