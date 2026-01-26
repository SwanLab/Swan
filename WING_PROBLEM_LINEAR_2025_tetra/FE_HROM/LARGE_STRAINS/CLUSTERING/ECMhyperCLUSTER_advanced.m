function [ECMdata_cluster,setCandidates,maX_INTERNAL_FORCES] = ECMhyperCLUSTER_advanced(BasisU_cluster,BasisPone_cluster,BasisStwo_cluster,...
    OPERFE,DISP_CONDITIONS,DATA,DATAoffline)



DATAoffline = DefaultField(DATAoffline,'TryUseSameECMpointsForAllClusters',1) ; % 

if DATAoffline.TryUseSameECMpointsForAllClusters ==1 
    % Sequential ECM 
    [ECMdata_cluster,setCandidates,maX_INTERNAL_FORCES] = ECMhyperCLUSTER_sequential(BasisU_cluster,BasisPone_cluster,BasisStwo_cluster,...
    OPERFE,DISP_CONDITIONS,DATA,DATAoffline) ; 
    
else
    error('Option not implemented')
    
end 

