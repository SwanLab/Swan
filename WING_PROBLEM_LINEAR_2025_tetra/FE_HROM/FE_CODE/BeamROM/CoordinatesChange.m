function NewCOOR = CoordinatesChange(OldCoor,OldCenter,NewCenter,RotationMatrix)
% 
% Change of coordinates
% % JAHO, 31-Dec-2017
%

NewCOOR  = zeros(size(OldCoor)) ;
for idim=1:size(OldCoor,1) % Relative coordinates (to the OldCenteroid)
    NewCOOR(idim,:) = OldCoor(idim,:) - OldCenter(idim) ;
end
% Rotated coordinates
NewCOOR = RotationMatrix*NewCOOR ;
% Translation
for idim = 1:size(OldCoor,1)
    NewCOOR(idim,:) = NewCOOR(idim,:) + NewCenter(idim) ;
end