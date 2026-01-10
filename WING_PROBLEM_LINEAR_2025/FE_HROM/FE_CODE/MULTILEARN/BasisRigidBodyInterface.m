function [BasisINTrb COORrefPOINT]=  BasisRigidBodyInterface(COORref,NODESfaces,DATAINM)

COOR_FACE = COORref(NODESfaces{1},:) ;
    COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
    COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
    BasisINTrb = ConstructBasisRigidBody(COORrel,DATAINM) ; %