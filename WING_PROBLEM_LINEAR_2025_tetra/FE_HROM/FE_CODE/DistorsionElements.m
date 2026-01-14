function DATAOUT = DistorsionElements(COOR,CN,TypeElement) 
% Function for computing the distorsion of FE elements 
% JAHO, 3-May-2020, 50th of confinment (COVID-19)
% --------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
nnodeE = size(CN,2) ; 
switch TypeElement 
    case 'Hexahedra' 
        
            LINES_ANGLES = {[1,2],[2,3],[3,4],[4,1],...
                            [5,6],[6,7],[7,8],[8,5],...
                            [1,5],[2,6],[3,7],[4,8]} ; 
             CNsides = cell2mat(LINES_ANGLES');        
             REF_angle = 90 ; 
              PERMT_ANGLES = {[1,2],[1,3],[2,3]} ;
    case 'Quadrilateral'
         LINES_ANGLES = {[1,2],[2,3],[3,4],[4,1]} ; 
             CNsides = cell2mat(LINES_ANGLES');        
             REF_angle = 90 ; 
              PERMT_ANGLES = {[1,2]} ;
    
    otherwise
       %  error('Not implemented')
       LINES_ANGLES = [] ; 
end

if isempty(LINES_ANGLES) 
  DATAOUT = [] ,
    
else

ndim = size(COOR,2) ; 
nelem = size(CN,1) ; 
maxNODES = max(CNsides(:)) ;
nangles = maxNODES*ndim ; 
ANGLES_avg = zeros(nelem,1) ; 
ANGLE_max =  zeros(nelem,1) ; 
ANGLE_min = 1e20*ones(nelem,1) ; 
distorsionANGLE = zeros(nelem,1) ;  
for inode1 = 1:maxNODES
    [elemLOC,nodeLOC] =  find(CNsides == inode1) ;
    inode1_GLO = CN(:,inode1) ;    
    X1 = COOR(inode1_GLO,:) ;
    vLINE = cell(1,length(elemLOC)) ; 
    for kkk = 1:length(elemLOC)
        icmpl =   setdiff(1:2,nodeLOC(kkk)) ;
        inode2 = CNsides(elemLOC(kkk),icmpl) ;   
        inode2_GLO = CN(:,inode2) ;    
        % Coordinates node
        X2 = COOR(inode2_GLO,:) ;
        dx = X1-X2;
        ndx = sqrt(sum(dx.^2,2)) ;
        vLINE{kkk} = zeros(nelem,ndim) ;  
        for idim = 1:ndim
        vLINE{kkk}(:,idim) = dx(:,idim)./ndx ;    
        end
    end
    
   
    for iii = 1:length(PERMT_ANGLES)
        combLINES = PERMT_ANGLES{iii}; 
        vect1=  vLINE{combLINES(1)} ; 
        vect2 =vLINE{combLINES(2)} ; 
        pESC = zeros(nelem,1) ; 
        
        for idim = 1:ndim 
            pESC = pESC + vect1(:,idim).*vect2(:,idim) ; 
        end
  
        newANGLE = acosd(pESC) ; 
        ANGLES_avg = ANGLES_avg +newANGLE ; 
        ANGLE_min = min([ANGLE_min,newANGLE],[],2) ; 
          ANGLE_max = max([ANGLE_max,newANGLE],[],2) ; 
          distorsionANGLE = distorsionANGLE + abs(newANGLE-REF_angle) ; 
    end
    
    
end

DATAOUT.ANGLES_avg = ANGLES_avg/nangles ; 
DATAOUT.distorsionANGLE = distorsionANGLE/nangles ; 
DATAOUT.ANGLE_max = ANGLE_max; 
DATAOUT.ANGLE_min = ANGLE_min; 

end
