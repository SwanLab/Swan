function ModesMatrix = CollectModesSeveralProjects_PCNO(NAME_MODES,FOLDER,TYPE,DATAINPUT)

ModesMatrix = [] ;


for ifile = 1:length(NAME_MODES)
    load([FOLDER,NAME_MODES{ifile}],'BASES','DATA_REFMESH') ;
    % COOR{ifile} = DATA_REFMESH.COOR ;
    %     xMAX = max(COOR{ifile}(:,1)) ;
    %     xMIN = min(COOR{ifile}(:,1)) ;
    %     LENGTHS(ifile) = xMAX-xMIN ;
    Basis = BASES.(TYPE).U ;
    
   
  
    
      if DATAINPUT.NORMALIZE_MODES_BEFORE_SVD == 1
          
          for iii = 1:size(Basis,2)
          normBasis =norm(Basis(:,iii)) ; 
          Basis(:,iii) = Basis(:,iii)/normBasis  ; 
          end
          
      end
    
      if DATAINPUT.ESTIMATE_ALL_MODES_AT_ONCE  ==1 
             Basis = Basis(:) ; 
        
    else
         Basis = Basis(:,DATAINPUT.imode) ;
    end
    

    
    
    
    ModesMatrix = [ModesMatrix,Basis] ;
    
end