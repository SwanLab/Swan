function BasisRrb = ConstructBasisRigidBody(COORrel,DATAIN)

if nargin == 1
    DATAIN.ORTHOGONAL_RIGID_BODY_MODES = 0 ;
end

ndim = size(COORrel,2) ;

if ndim == 2
    BasisRrb =  zeros(prod(size(COORrel)),3) ;
    BasisRrb(1:ndim:end,1) = 1;
    BasisRrb(2:ndim:end,2) = 1;
    BasisRrb(1:ndim:end,3) = -COORrel(:,2);
    BasisRrb(2:ndim:end,3) = +COORrel(:,1);  % or the other way around ?
else
    
    BasisRrb =  zeros(prod(size(COORrel)),6) ;
    BasisRrb(1:3:end,1) = 1 ;
    BasisRrb(2:3:end,2) = 1;
    BasisRrb(3:3:end,3) = 1;
    
    % ------------------------------------------------
    % CHANGE INTRODUDUCED 3th-october-2021
    %------------------------------------
    % OLD VERSION 
    % ---------------------------
%      % Rotation modes
%     BasisRrb(2:3:end,4) = COORrel(:,3);
%     BasisRrb(3:3:end,4) = -COORrel(:,2);
%     
%     BasisRrb(1:3:end,5) = -COORrel(:,3);
%     BasisRrb(3:3:end,5) = COORrel(:,1);
%     
%     BasisRrb(1:3:end,6) = COORrel(:,2);
%     BasisRrb(2:3:end,6) = -COORrel(:,1);


% NEW VERSION 
% ----------------------------
% See  /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/MultiECMgen/DiscussionRIGIDBODY.xoj  /*.tex 
    % Rotation modes
    % [R_r]_I   = ROtation modes of node I, of coordinates [X]_I = [x1,x2,x3]
    %                MODE  4   5   6
    %                ------------------
    % = -spin([X]_I)  =  [0  +x3 -x2 
    %                      -x3  0  +x1
    %                      +x2  -x1  0] = 
    
    BasisRrb(2:3:end,4) = -COORrel(:,3);
    BasisRrb(3:3:end,4) = +COORrel(:,2);
    
    BasisRrb(1:3:end,5) = +COORrel(:,3);
    BasisRrb(3:3:end,5) = -COORrel(:,1);
 
     BasisRrb(1:3:end,6) = -COORrel(:,2);
     BasisRrb(2:3:end,6) = +COORrel(:,1);
   

% ---------------------------------------------------------------------








end

DATAIN = DefaultField(DATAIN,'ORTHOGONAL_RIGID_BODY_MODES',0) ;

if DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 1
    BasisRrb = orth(BasisRrb)  ;
elseif DATAIN.ORTHOGONAL_RIGID_BODY_MODES== 2
    NNN = sqrt(sum(BasisRrb.^2,1)) ;
    for imode = 1:size(BasisRrb,2)
        BasisRrb(:,imode) = BasisRrb(:,imode)/NNN(imode) ;
    end
end