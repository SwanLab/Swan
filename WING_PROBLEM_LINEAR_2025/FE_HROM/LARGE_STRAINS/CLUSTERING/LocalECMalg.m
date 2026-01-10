function [ECMdata,setCandidates ]= LocalECMalg(wSTs,DATAoffline,Q,DATA,setCandidates)

if nargin == 0
    load('tmp.mat')
end
 
% -------------------------
DATA_ECM = [] ;
DATAoffline = DefaultField(DATAoffline,'errorECM',0) ; 
DATA_ECM.TOL = DATAoffline.errorECM ;
if isempty(setCandidates)
    DATA_ECM.IND_POINTS_CANDIDATES = 1:length(wSTs) ;
else
    DATA_ECM.IND_POINTS_CANDIDATES = setCandidates ;
end
[setPoints,wRED,ERROR_GLO,DATAOUT]= EmpiricalCubatureMethod_CANDcompl(Q,[],wSTs,DATA_ECM)  ;

setCandidates = unique([setCandidates;setPoints]) ;
 
ECMdata.setPoints = setPoints ;
ECMdata.wRED = wRED ;
proporP = length(setPoints)/length(wSTs)*100;
disp(['Number of ECM points = ',num2str(length(setPoints)),' (',num2str(proporP),' % total)'])
disp(['List of selected ECM points  = ',num2str(length(setPoints))])

 