function angCOV = CorrelationAnglesPairFaces(Drb) 



iM =  1:size(Drb,1)/2 ;  % Entries corresponding to sides  1-2
iP =  (size(Drb,1))/2+1:size(Drb,1) ;  % Entries corresponding to sides 3-4

DrbP = Drb(iM,:) ; % 
DrbM = Drb(iP,:) ;

% Degree of linear correlation (ANGLE) between pair of modes 
% ----------------------------------------------------------------------
for imode = 1:size(DrbP,1) 
    DrbP(imode,:) = DrbP(imode,:)/norm(DrbP(imode,:)) ; 
    DrbM(imode,:) = DrbM(imode,:)/norm(DrbM(imode,:)) ; 
end 
% COVARIANCE MATRIX RIGID BODY MODES, POSITIVE SIDES, AND NEGATIVE SIDES (1-2 / 3-4)
% -----------------
COV = DrbP*DrbM' ; 
COV = diag(COV) ; 
% Angle formed by pair of modes 
angCOV = acos(abs(COV))*180/pi; 