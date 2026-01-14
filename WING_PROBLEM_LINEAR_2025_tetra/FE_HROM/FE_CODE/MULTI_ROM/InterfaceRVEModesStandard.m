function  [DATAOUT,V,TEXTP] =  InterfaceRVEModesStandard(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
    SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS)


DATAIN = DefaultField(DATAIN,'INCLUDE_FLUCTUATION_MODES_FACES',1) ;
% If = 0, fluctuation modes are ignored and
% interface modes are rigid
% body modes
if  DATAIN.INCLUDE_FLUCTUATION_MODES_FACES==0
    % Rigid body modes for the interfaces
    % ------------------------------------
    DATAOUT.BasisInt = Vrb ;
    V = Vrb ;
else
    RotationMatrixALL = cell(1,nfaces) ;

    V = cell(size(Vrb)) ;
    %   RotatedReactions = cell(size(Vrb)) ;
    %   RotatedReactionsT = cell(1,length(FACES_GROUPS)) ;
    for ifgroup=1:length(FACES_GROUPS)
        f1 = fI{FACES_GROUPS{ifgroup}(1)} ;
        f2 = fI{FACES_GROUPS{ifgroup}(2)} ;
        iface = FACES_GROUPS{ifgroup}(1) ;
        [Vall,RotationMatrixLOC,TEXTP ] = ...
            ReactionAndInterfaceLocalModes_RVE(BasisUdef,BasisRdef,f1,f2,...
            DATAIN.TOL_SINGULAR_VALUES_Hqr,...
            Vrb{iface},M{iface},DATAIN,SinvVal_Udef,SinvVal_Rdef,ifgroup,TEXTP)  ;
        V{FACES_GROUPS{ifgroup}(1)} = Vall ;
        V{FACES_GROUPS{ifgroup}(2)} = Vall ;
        if ~isempty(RotationMatrixLOC)
            RotationMatrixALL{FACES_GROUPS{ifgroup}(1) } = RotationMatrixLOC{1} ;
            RotationMatrixALL{FACES_GROUPS{ifgroup}(2) }  = RotationMatrixLOC{2} ;
        end        
    end   
    DATAOUT.BasisInt = V  ; 
end