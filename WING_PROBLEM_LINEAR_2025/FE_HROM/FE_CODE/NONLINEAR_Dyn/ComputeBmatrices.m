function [wSTs,  XeALL,  Bst, BstW, wST,  DATA] = ...
    ComputeBmatrices(COOR,CN,TypeElement,DATA)
if nargin == 0
    load('tmp0.mat')
end

disp(['DATA.RECALCULATE_STIFFNESS=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
if isfield(DATA,'AREA')
    AREA = DATA.AREA ;
else
    AREA = [] ;
end
%% RENUMERING CONNECTIVITY MATRIX (to ensure small bandwidth of both BB and BBnw)


%   dbstop('53')

% Vectorized code

switch TypeElement
    
    case 'Quadrilateral'
        
        DATA = DefaultField(DATA,'USE_BBAR_ELEMENT',1) ;
        
        if  DATA.USE_BBAR_ELEMENT == 0 || size(CN,2) ~=4
            
            % GEneral method, standard elements. Faster method, every
            % operation is vectorized
            [ wSTs  XeALL  Bst wST DATA.posgp]=...
                ComputeBmatrices_loc(COOR,CN,TypeElement,DATA) ;
            
            if DATA.ComputeBstW == 1
            wSTdiag = diag(sparse(wST)) ;
            
            BstW = wSTdiag*Bst ;
            else
               BstW = [];  
            end
            
        else
            % BBar method, 2D. Not so fast ...
            [Bst,BstW,wSTs,wST,DATA.posgp] = BMatQuad_JAHO(COOR,CN,1,DATA);
            % COORDINATE MATRIX arranged in a nelem*ndim x nnodeE matrix
            XeALL= COORallELEM(size(COOR,2),size(CN,1),size(CN,2),CN,COOR) ;
        end
        
    otherwise
        
        % GEneral method, standard elements. Faster method, every
        % operation is vectorized
        [ wSTs  XeALL  Bst wST DATA.posgp]=...
            ComputeBmatrices_loc(COOR,CN,TypeElement,DATA) ;
        
        
        
        if  DATA.ComputeBstW == 1
        wSTdiag = diag(sparse(wST)) ;
        
        BstW = wSTdiag*Bst ;
        else
            BstW = [] ; 
        end
        
        
end



%    dbstop('61')
disp('Storing B Matrices...')
%     save(DATA.nameWORKSPACE,'K','CN','wSTs','XeALL','Cglo','Bst','wST','MaterialType','COOR','-append');

if isfield(DATA,'densGLO')
    density = DATA.densGLO ;
else
    density = [] ;
end


save(DATA.nameWORKSPACE,'wST','wSTs','XeALL',...
    'AREA','density','Bst','BstW','-append');


