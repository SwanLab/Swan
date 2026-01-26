function  ECMdata = SAW_ECM_large1param(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT,DATA_interp)
 
if nargin == 0
    load('tmp3.mat')
    DATAoffline.ECM_Number_Snapshots_Overlapping =0; 
%     DATAoffline.SetCandidates_given_by_the_userECM = 72;  
    DATAoffline.IntegrateExactlyVolume_SAW_ECM  = 1 ; 
    DATAoffline.errorFINT = 1e-4; 
  close all
end


DATAoffline = DefaultField(DATAoffline,'ECM_Number_Snapshots_Overlapping',0) ; % =0; 


if DATAoffline.ECM_Number_Snapshots_Overlapping >0 
  
 ECMdata = SAW_ECM_large1param_OVERLAP(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT) ;
 error('Option not maintained in the online phase')
else
    ECMdata = SAW_ECM_large1param_CONT(DATAoffline,SNAPfint,DATA,OPERFE,qLATENT,DATA_interp) ;
end