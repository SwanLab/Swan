function COOR = ScalingSliceLoc(COOR,Cf1,SCALE_FACTOR)

% ----------------------------------
COORrel = zeros(size(COOR)) ; 
% Relative coordinates
for idim = 1:length(Cf1)
    COORrel(:,idim) =  COOR(:,idim)- Cf1(idim) ; 
end
% Scaling 
COOR = zeros(size(COORrel)) ; 
FIELDS = {'X','Y','Z'};
for idim = 1:length(Cf1)
    e = SCALE_FACTOR.(FIELDS{idim}) ; 
    COOR(:,idim) = Cf1(idim) + COORrel(:,idim)*e ; 
end