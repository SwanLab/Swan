function    [RESULTANTS,DISP_x]= ReactionsResultantDISP(SNAPreactions,DOFr,MESH,DATAINPUTloc,Fbody_U,DATA,DISP_LOC_HROM,...
        BasisUall,Finer, Fdamp, FintR)

if nargin == 0
    load('tmp3.mat')
end


ndim = size(MESH.COOR,2) ;
WEIGHT  =  norm(sum(reshape(Fbody_U,ndim,[]),2)) ;
disp(['Total weight =',num2str(WEIGHT),' s.u.'])


DATAINPUTloc = DefaultField(DATAINPUTloc,'SURFACE_TO_STUDY_REACTIONS',1) ;
DATAINPUTloc = DefaultField(DATAINPUTloc,'SURFACE_TO_STUDY_DISPLACEMENTS',2) ;


[DOFr,aaa ]= sort(DOFr) ;
SNAPreactions = SNAPreactions(aaa,:);


isurf = DATAINPUTloc.SURFACE_TO_STUDY_REACTIONS ;
isurfDISP = DATAINPUTloc.SURFACE_TO_STUDY_DISPLACEMENTS ;

% Centroid of these surface, as well as associated nodes
% ------------------------------------------------------
NODES = sort(MESH.NODES_FACES{isurf}) ;
NODES_DISP = sort(MESH.NODES_FACES{isurfDISP}) ;

CENTROID = MESH.PROPERTIES_FACES{isurf}.CENTROID' ;


% Coordinates mass center with respect to rotation center
CENTER_mass_rel = MESH.CENTER_MASS' -CENTROID ;

% Maximum moment is achieved when the position vector of the center of mass
% relative to the rotation center and gravity load are normal
Moment_mass = norm(CENTER_mass_rel)*WEIGHT ;


% Associated DOFs
DOFr_surf = small2large(NODES,ndim) ;
DOFr_surf_DISP = small2large(NODES_DISP,ndim) ;


[DOFr_intersect,iALL,iLOC] = intersect(DOFr,DOFr_surf) ;
if length(DOFr_intersect) ~= length(DOFr_surf)
    error('Incomplete set of DOFs')
end


[DOFr_intersect_DISP,iALL_DISP,iLOC] = intersect(DOFr,DOFr_surf_DISP) ;
if length(DOFr_intersect) ~= length(DOFr_surf)
    error('Incomplete set of DOFs')
end

% Check that the associated DOFr are complete, in the sense that
REACTS = SNAPreactions(iALL,:) ;
Finer = Finer(iALL,:) ;

Finer = sum(Finer(1:2:end,:),1) ; 
 
Fdamp = Fdamp(iALL,:) ;
Fdamp = sum(Fdamp(1:2:end,:),1) ; 


FintR = FintR(iALL,:) ;
FintR = sum(FintR(1:2:end,:),1) ; 

Ftotal = sum(REACTS(1:2:end,:),1) ; 

% figure(788)
% hold on 
% h1 = plot(Finer)
% h2 = plot(Fdamp)
% h3 = plot(FintR)
% h4 = plot(Ftotal)
% 
% legend([h1,h2,h3,h4],{'Finer','Fdamp','FintR','Ftotal'})


 

DISP = BasisUall(DOFr_surf_DISP,:)*DISP_LOC_HROM ;  
DISP_x = sum(DISP(1:2:end,:),1) ;
DISP_x = DISP_x/length(DOFr_surf_DISP)*2 ; 

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
