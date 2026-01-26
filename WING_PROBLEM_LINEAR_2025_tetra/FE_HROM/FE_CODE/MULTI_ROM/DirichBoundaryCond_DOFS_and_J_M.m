
function [J,DOFS,LENGTH,M] = DirichBoundaryCond_DOFS_and_J_M(MESH2D,CNb,ielem,Js1,Js2,ndim,xC,yC)

if nargin == 0
    load('tmp.mat')
end

NodMid= intersect(MESH2D.NODESmid,CNb(ielem,:));
CNbVERT= setdiff(CNb(ielem,:),MESH2D.NODESmid ) ;
% Correspondence in local numbering
NODE(1) =  find(CNbVERT(1)==MESH2D.NODESvert) ;
NODE(2) =  find(CNbVERT(2)==MESH2D.NODESvert) ;

LENGTH = norm(MESH2D.COORvert(NODE(1),:)-MESH2D.COORvert(NODE(2),:));
COORmid = 0.5*(MESH2D.COORvert(NODE(1),:)+MESH2D.COORvert(NODE(2),:)) ; 

%% Quadrilateral element containing NODES1 and NODE2
OK_1 = zeros(size(MESH2D.CN,1)) ;
OK_2 = OK_1 ;
for  jnode = 1:size(MESH2D.CN,2)
    OK_1(:,jnode) =  (NODE(1) ==MESH2D.CN(:,jnode));
    OK_2(:,jnode) =  (NODE(2) ==MESH2D.CN(:,jnode));
end

OK = sum(OK_1,2).*sum(OK_2,2) ;
INDEX_ELEM = find(OK==1) ;
INX_LOC_1 = find(OK_1(INDEX_ELEM,:)==1) ;
INX_LOC_2 = find(OK_2(INDEX_ELEM,:)==1) ;

%% Type of element
itype = MESH2D.MaterialType(INDEX_ELEM) ;
J1 = Js1{itype} ;
J2 = Js2{itype} ;

FACES = {[1 2],[2 3 ],[4 3],[1,4]} ;
% Indefication of faces
iface= 1 ;

while iface <=  length(FACES)
    CONJD = setdiff(FACES{iface},[INX_LOC_1 INX_LOC_2]) ;
    if isempty(CONJD)
        
        break
    end
    iface= iface + 1;
end

if iface == 1 | iface == 3
    J = J1 ;
        % Matrix relating input data displacements with generalized
        % displacement at the studied line
        nmodes = 6  ;
        M = diag(ones(nmodes,1)) ;
        % Midside node, coordinates
        y_0  = COORmid(2)-yC ;
        M(3,4) = y_0 ;
        M(1,6) = -y_0 ;
    
 elseif iface == 2  | iface == 4
    J = J2 ;
        % Matrix relating input data displacements with generalized
    % displacement at the studied line  
    nmodes = 6  ; 
    M = diag(ones(nmodes,1)) ; 
    % Midside node, coordinates 
    x_0  = COORmid(1)-xC ; 
    M(2,6) = x_0 ; 
    M(3,5) = -x_0 ;  
else
    error('Check the implementation....')
end

NODE_ORD = NODE;
if INX_LOC_1 ~= FACES{iface}(1)
    NODE_ORD(1) = NODE(2) ;
    NODE_ORD(2) = NODE(1) ;
end
DOFS = small2large(NODE_ORD,ndim) ;

end


