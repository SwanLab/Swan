function ModesMatrix = CollectModesProjects(NAME_MODES,FOLDER,TYPE,DATAINPUT)

if nargin == 0
    load('tmp1.mat')
end 
ModesMatrix = {} ;

DATAINPUT = DefaultField(DATAINPUT,'NO_PHYSICAL_CORRESPONDENCE_BETWEEN_MODES',0) ;

for ifile = 1:length(NAME_MODES)
    load([FOLDER,NAME_MODES{ifile}],'BASES','DATA_REFMESH') ;
    % COOR{ifile} = DATA_REFMESH.COOR ;
    %     xMAX = max(COOR{ifile}(:,1)) ;
    %     xMIN = min(COOR{ifile}(:,1)) ;
    %     LENGTHS(ifile) = xMAX-xMIN ;
    switch  TYPE
        case {'BRED'}
            load([FOLDER,NAME_MODES{ifile}],'BASES','Bdom') ;
            Udisp =  BASES.DISPLACEMENTS.U ; 
             ModesMatrix{ifile} = Bdom*Udisp ;
        otherwise
            ModesMatrix{ifile} = BASES.(TYPE).U ;
    end
    
%     if ifile == 1
%         BasisREF = Basis ; % Set of reference modes
%     else
%         if    DATAINPUT.NO_PHYSICAL_CORRESPONDENCE_BETWEEN_MODES == 1
%             [Basis,~, ~] = SVDT(Basis) ;
%             if DATAINPUT.ESTIMATE_ALL_MODES_AT_ONCE  ==1
%                 error('Non-compatible options')
%             end
%         end
%     end
%     if DATAINPUT.NORMALIZE_MODES_BEFORE_SVD == 1
%         
%         for iii = 1:size(Basis,2)
%             normBasis =norm(Basis(:,iii)) ;
%             Basis(:,iii) = Basis(:,iii)/normBasis  ;
%         end
%         
%     end
    
%     if DATAINPUT.ESTIMATE_ALL_MODES_AT_ONCE  ==1
%         Basis = Basis(:) ;
%     else
%         
%         
%         if DATAINPUT.NO_PHYSICAL_CORRESPONDENCE_BETWEEN_MODES == 1
%             Basis = Basis*(Basis\BasisREF(:,DATAINPUT.imode)) ;
%             Basis = Basis/norm(Basis) ;
%         else
%             % Standard approach
%             Basis = Basis(:,DATAINPUT.imode) ;
%         end
%         
%         
%     end
%     
    
    
%     
%     
%     ModesMatrix = [ModesMatrix,Basis] ;
    
end