function    [RESULTANTS,ANG,WEIGHT,Moment_mass ]= ReactionsResultant(SNAPreactions,DOFr,MESH,DATAINPUTloc,Fbody_U,DATA)

if nargin == 0
    load('tmp3.mat')
end


ndim = size(MESH.COOR,2) ;
WEIGHT  =  norm(sum(reshape(Fbody_U,ndim,[]),2)) ;
disp(['Total weight =',num2str(WEIGHT),' s.u.'])


DATAINPUTloc = DefaultField(DATAINPUTloc,'SURFACE_TO_STUDY_REACTIONS',1) ;

[DOFr,aaa ]= sort(DOFr) ;
SNAPreactions = SNAPreactions(aaa,:);


isurf = DATAINPUTloc.SURFACE_TO_STUDY_REACTIONS ;

% Centroid of these surface, as well as associated nodes
% ------------------------------------------------------
NODES = sort(MESH.NODES_FACES{isurf}) ;
CENTROID = MESH.PROPERTIES_FACES{isurf}.CENTROID' ;


% Coordinates mass center with respect to rotation center
CENTER_mass_rel = MESH.CENTER_MASS' -CENTROID ;

% Maximum moment is achieved when the position vector of the center of mass
% relative to the rotation center and gravity load are normal
Moment_mass = norm(CENTER_mass_rel)*WEIGHT ;


% Associated DOFs
DOFr_surf = small2large(NODES,ndim) ;

[DOFr_intersect,iALL,iLOC] = intersect(DOFr,DOFr_surf) ;
if length(DOFr_intersect) ~= length(DOFr_surf)
    error('Incomplete set of DOFs')
end

% Check that the associated DOFr are complete, in the sense that
REACTS = SNAPreactions(iALL,:) ;
%REACTS = reshape(REACTS,ndim,[]) ;
nsteps = size(SNAPreactions,2) ;

COORr = MESH.COOR(NODES,:)' ;
COORrREL = zeros(size(COORr)) ;
for idim=1:ndim
    COORrREL(idim,:) = COORr(idim,:)-CENTROID(idim) ;
end


if ndim == 2
    RESULTANTS = zeros(3,nsteps) ;
    COORrREL = [COORrREL; zeros(1,size(COORrREL,2))]  ;
else
    RESULTANTS = zeros(6,nsteps) ;
end


DATAINPUTloc = DefaultField(DATAINPUTloc,'AngleVersusTime_vector',[]) ;

for istep = 1:nsteps
    % This may be optimized...
    % ------------------------
    %     if istep == 1020
    %         disp('borrar esto')
    %     end
    REACTSloc = reshape(REACTS(:,istep),ndim,[]) ;
    
    
 %   if  isempty(DATAINPUTloc.AngleVersusTime_vector)
        t = DATA.STEPS(istep) ;
  %      ANG = DATAINPUTloc.ROTATION_ANGLE_time(t) ;
        
   % else
    %    ANG = DATAINPUTloc.AngleVersusTime_vector(istep) ;
    %end
    
    
    
    if ndim == 2
    %    R= [cosd(ANG) -sind(ANG)
     %       sind(ANG)  cosd(ANG)] ;
     %   REACTSloc_ROT = R'*REACTSloc ;
        RESULTANTS(1:2,istep) = sum(REACTSloc,2) ;
        
        Moments = cross(COORrREL,[REACTSloc; zeros(1,size(REACTSloc,2))]) ;
        RESULTANTS(3,istep) = sum(Moments(3,:),2) ;
        
    else
        
        
%         R= [cosd(ANG) -sind(ANG)  0
%             sind(ANG)  cosd(ANG)  0
%             0           0        1] ;
     %   REACTSloc_ROT = R'*REACTSloc ;
        RESULTANTS(1:3,istep) = sum(REACTSloc,2) ;
        
        Moments = cross(COORrREL,REACTSloc) ;
        RESULTANTS(4:6,istep) = sum(Moments(1:3,:),2) ;
        
        
    end
    
end

% if isempty(DATAINPUTloc.AngleVersusTime_vector)
%     ANG = DATAINPUTloc.ROTATION_ANGLE_time(DATAINPUTloc.TIME_STEPS) ;
% else
%     ANG = DATAINPUTloc.AngleVersusTime_vector ;
% end


%





%
% DOFrORDER = sort(DOFr) ;
%
%
% ReactSORTED = React(DOFrORDER);
% ReactSORTED = reshape(ReactSORTED,ndim,[]) ;
%
% DIV = DOFr/ndim ;
% ListNodes = ceil(DIV) ;
% ListNodes = unique(ListNodes) ;
%
% COORr = COOR(ListNodes,:)' ;
% COORrREL = zeros(size(COORr)) ;
% for idim=1:ndim
%     COORrREL(idim,:) = COORr(idim,:)-CENTROID(idim) ;
% end
%
% Moments = cross(COORrREL,ReactSORTED) ;
%
% MomentsResult = sum(Moments,2) ;
%
%
% Resultants = sum(ReactSORTED,2) ;
