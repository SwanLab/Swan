function [NORMALSv,TANGENTSv]= NormalsBoundaryLocal(COOR,CNb)
% Compute normals of list CNb (connectivities)

if nargin == 0
    load('tmp.mat')
end

if isempty(CNb)
    error('You must generathe the Boundary Mesh !!!')
end

ndim = size(COOR,2) ;
NORMALSv = zeros(ndim,size(CNb,1));

Xloc1 = COOR(CNb(:,1),:)';
Xloc2 = COOR(CNb(:,2),:)';
t1 = Xloc2-Xloc1 ; % Tangential vector

normt1 = sqrt(sum(t1.^2,1)) ;
for idim = 1:ndim
    t1(idim,:)  = t1(idim,:)./normt1 ;
end


if ndim == 3
    
    Xloc3 = COOR(CNb(:,3),:)';
    t2 = Xloc3-Xloc1 ; % Tangential vector
    
    normt2 = sqrt(sum(t2.^2,1)) ;
    for idim = 1:ndim
        t2(idim,:)  = t2(idim,:)./normt2 ;
    end
    
    n = cross(t1,t2) ;
    
else
    n = zeros(size(t1)) ;
    n(1,:) = -t1(2,:) ;
    n(2,:) = t1(1,:) ;
    
end

normN = sqrt(sum(n.^2,1)) ;
for idim = 1:size(n,1)
    NORMALSv(idim,:)  = n(idim,:)./normN ;
end

% Now we have to enforce that all normals point to the same DIRECTION (see KW:NormalsSameDirections)


% REFERENCE VECTOR
NORMALSv_ref = NORMALSv(:,1) ;
% Scalar product (they should be all positive, = 1 if its a plane )
PROD_ESC = zeros(1,size(NORMALSv,2)) ;
%
for idim = 1:ndim
    PROD_ESC = PROD_ESC + NORMALSv_ref(idim)*NORMALSv(idim,:) ;
end
%
III = find(PROD_ESC <1) ;

if ~isempty(III)
    % Two things might be happening: 1) Curved surface, embracing more
    % than 90 degrees (see KW:NormalsSameDirections)
    % 2) Abrupt change due to a wrong choice of normals
    
    % Assuming is the second case, we proceed by changing the orientation
    % of the normal
    for idim = 1:ndim
        NORMALSv(idim,III) = -NORMALSv(idim,III) ;
    end
    
end

% TANGENTIAL VECTORS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEe DOCS, kw:TANGENTIALvectors
if ndim == 2
    
    TANGENTSv = [NORMALSv(2,:); - NORMALSv(1,:)];
    
    
else
    
    % Check if the tangent plane is either XY,XZ or YZ
    FOUND = 0 ;
    while idim <=ndim
        TOL = 1e-8 ;
        CHECK = abs(NORMALSv(idim,:))<TOL ;
        if all(CHECK)
            % A tangent vector is
            tangent1 = zeros(ndim,1) ;
            tangent1(idim) = 1;
            
            % Therefore
            tangent1 = repmat(tangent1,1,size(NORMALSv,2)) ;
            tangent2 = cross(NORMALSv,tangent1) ;
            FOUND=1;
            break
        end
        idim = idim + 1;
    end
    
    if FOUND ==0
        % The tangent plane is none of the planes XY, YZ or XZ
        %-----------------------------------------------------
        % We take a reference element, determine arbitrarily its tangent
        % vectors, and then choose the tangent vectors of the remaining
        % elements so that the alignment with respect to the reference one
        % is maximized (projection )
        iref  = 1;
        t1REF = t1(:,iref) ;
        % Now we have to find the projection of t1REF in the plane spanned
        % by the remaining t1 and t2.  KW:MaximimizeAlignmentTangent
        
        C1 =  zeros(1,size(t1,2)) ;
        for  idim = 1:ndim
            C1 = C1 + t1(idim,:).*t1REF(idim) ;
        end
        U1 = zeros(size(t1)) ;
        for idim = 1:ndim
            U1(idim,:) = t1(idim,:).*C1  ;
        end
        
        C2 =  zeros(1,size(t2,2)) ;
        for  idim = 1:ndim
            C2 = C2 + t2(idim,:).*t1REF(idim) ;
        end
        U2 = zeros(size(t2)) ;
        for idim = 1:ndim
            U2(idim,:) = t2(idim,:).*C1  ;
        end
        
        tangent1 = U1+U2;
        
        
        tangent2 = cross(NORMALSv,tangent1) ;
        
        
    end
    
    
    
    TANGENTSv{1} = tangent1 ;
    TANGENTSv{2} = tangent2 ;
    
    
end







PRINT_LOCAL = 0;

if PRINT_LOCAL ==1
    figure(590)
    hold on
    axis equal
    plot(COOR(:,1),COOR(:,2),'r*') ;
    tiN = norm(t1(:,1)) ;
    FACTOR = 3;
    NORMALSv = FACTOR*NORMALSv*tiN ;
    
    for ielem = 1:size(CNb,1)
        DIFFCOOR = 0.5*(COOR(CNb(ielem,2),:) + COOR(CNb(ielem,1),:)) ;
        
        xINI = DIFFCOOR' ;
        xFIN = DIFFCOOR' + NORMALSv(:,ielem) ;
        
        plot([xINI(1) xFIN(1)],[xINI(2), xFIN(2)],'r') ;
    end
    
end




end