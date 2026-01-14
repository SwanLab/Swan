function [K, CN, wSTs,  XeALL, Cglo, Bst, wST, MaterialType,  DATA] = ...
    StiffnessMatrixGlobal(COOR,CN,TypeElement,celasglo,DATA,MaterialType)

DATA  =DefaultField(DATA,'RECALCULATE_STIFFNESS',1) ;

disp(['DATA.=',num2str(DATA.RECALCULATE_STIFFNESS)]) ;
if isfield(DATA,'AREA')
    AREA = DATA.AREA ;
else
    AREA = [] ;
end
%% RENUMERING CONNECTIVITY MATRIX (to ensure small bandwidth of both BB and BBnw)
DATA = DefaultField(DATA,'RENUMBERED_OUTSIDE',0) ; 

DATA = DefaultField(DATA,'DO_NOT_STORE_STIFFNESS_AND_B_MATRIX',0) ; 
DATA = DefaultField(DATA,'VECTcode',1) ; 

if  DATA.RECALCULATE_STIFFNESS == 1
    %   dbstop('53')
    if DATA.VECTcode  == 0
        % Standard way (elementwise)
        K = ComputeK(COOR,CN,TypeElement, celasglo) ;
        wSTs= [] ; wsT = [] ; XeALL = [] ;Cglo = [] ; Bst = [] ;
        IndicesRenumberingElements = [] ;
    else
        % Vectorized code
        [K CN wSTs  XeALL Cglo Bst wST MaterialType IndicesRenumberingElements]=...
            ComputeKvect(COOR,CN,TypeElement, celasglo,DATA,MaterialType) ;
    end
    %    dbstop('61')
    DATA = DefaultField(DATA,'STORE_STIFFNESS',0) ; 
    if (DATA.STORE_STIFFNESS ==1 |  DATA.STORE_STIFFNESS ==2) %& DATA.DO_NOT_STORE_GEOMETRIC_INFORMATION ==0
        disp('Storing Stiffness Matrix...')
        %     save(DATA.nameWORKSPACE,'K','CN','wSTs','XeALL','Cglo','Bst','wST','MaterialType','COOR','-append');
        
        if isfield(DATA,'densGLO')
            density = DATA.densGLO ;
        else
            density = [] ;
        end
        
        
        save(DATA.nameWORKSPACE,'wST','wSTs','XeALL',...
            'AREA','density','-append');
        
        if DATA.RENUMBERED_OUTSIDE == 0
             save(DATA.nameWORKSPACE,'MaterialType','CN',...
            'IndicesRenumberingElements','-append');            
        end
        
        if DATA.DO_NOT_STORE_STIFFNESS_AND_B_MATRIX ==0 
            save(DATA.nameWORKSPACE,'Bst','K','Cglo','-append')
        end
        disp('Done')
        
        
    end
    
else
    
    DATA = DefaultField(DATA,'nameWORKSPACE_Kstiff',DATA.nameWORKSPACE) ;
    
    disp('Retrieving Stiffness Matrix...')
    %  dbstop('69')
    load(DATA.nameWORKSPACE_Kstiff,'K','CN','wSTs','XeALL','Cglo','Bst','wST','MaterialType');
    disp('Done')
end
