function  [U,S,V,h,h2,Ui_Si] = SVD_and_error(SNAP,nfigure,LEGENDG,NMODES_SHOW,COLOR,DATA,COLUMNS_RVEloc )

%dbstop('4')
if nargin ==  0
    load('tmp.mat')
    %  DATA.SVD_RANDOM = 1 ;
end
%SNAP = cell2mat(dRVE) ;
if iscell(SNAP)
    SNAP = cell2mat(SNAP) ;
end
h= [] ; h2 = [] ;

DATA = DefaultField(DATA,'SVD_RANDOM',1) ;
DATA = DefaultField(DATA,'TOL_LOC',0) ;

DATA = DefaultField(DATA,'REACTIONmodes_equal_DISPmodes',0) ;
DATA = DefaultField(DATA,'TOLERANCE_SVD_PROJECT',[] ) ; 

if ~isempty(DATA.TOLERANCE_SVD_PROJECT) &    (DATA.REACTIONmodes_equal_DISPmodes) ~= 0
    error('Incompatible options')

    
end
Ui_Si = [] ; 
DATA = DefaultField(DATA,'NMODES_PROJECT_LOC',[]) ; 
if   ~isempty(DATA.NMODES_PROJECT_LOC) 
     [nnnn mmmm]= cellfun(@size,COLUMNS_RVEloc) ;
     SNAP = mat2cell(SNAP,size(SNAP,1),mmmm);
    [U,S,V,Ui_Si] = RSVDcol(SNAP,DATA.NMODES_PROJECT_LOC) ; 
    eSVD = 0 ;

elseif (DATA.REACTIONmodes_equal_DISPmodes == 0 & isempty(DATA.TOLERANCE_SVD_PROJECT)) |  all(DATA.REACTIONmodes_equal_DISPmodes) 
    DATASVD.RELATIVE_SVD = 1 ;
    if DATA.SVD_RANDOM == 0
        [U,S,V,eSVD] = SVDT(SNAP,DATA.TOL_LOC,DATASVD) ;
    else
        [U,S,V,eSVD] = RSVDT(SNAP,DATA.TOL_LOC,[],0,DATASVD) ;
    end
    
elseif any(DATA.REACTIONmodes_equal_DISPmodes) 
   
        [U,S,V,eSVD] = SVD_selective(SNAP,DATA,COLUMNS_RVEloc) ; 

        
elseif ~isempty(DATA.TOLERANCE_SVD_PROJECT) 
    
    [U,S,V,eSVD] = SVD_accuracy_project(SNAP,DATA,COLUMNS_RVEloc) ; 
    

else
    error('Option not implemented')
    
end

DATA = DefaultField(DATA,'PLOT_RIGHT_SINGULAR_VECTORS',0); % .PLOT_RIGHT_SINGULAR_VECTORS =1 ; 
 
if DATA.PLOT_RIGHT_SINGULAR_VECTORS >0
    NMODES_RIGHT = min(DATA.PLOT_RIGHT_SINGULAR_VECTORS,length(S)) ;
    LLL = DATA.LEGEND_GRAPHS ;
    LLLeg = {} ;
    if DATA.LEGEND_GRAPHS(1) == 'D'
        nfigures = 400 ;
    else
        nfigures  = 401 ;
    end
    figure(nfigures)
    hold on
    xlabel('Domain')
    ylabel('Right S. Val.')
    colores = ColoresMatrix(NMODES_RIGHT );
    title(['Right singular v. x Sing. Val., ',DATA.LEGEND_GRAPHS])
    
    for iii = 1:NMODES_RIGHT
        hLOC(iii) = plot(V(:,iii)*S(iii),'Color',colores(iii,:)) ;
        LLLeg{iii} = ['i=',num2str(iii),',  s =',num2str(S(iii),4)] ; 
    end
    grid on
    legend(hLOC,LLLeg) ; 
    
    
end


if eSVD == 0  
    eSVD =  eps(1);
end

if ~isempty(S)  
    
    SingVsq =  (S.*S) ;
    nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
    SingVsq = sort(SingVsq);
    normEf2 = sqrt(cumsum(SingVsq)) ;
    normEf2 = sort(normEf2,'descend');
    % hold on
    % lS = log10(normEf2);
    % figure(1)
    % hold on
    % h = plot([1:length(S)-1],lS(2:end),'r');
    % xlabel('MODES')
    % ylabel('LOG10 SVD ERROR')
    % AAA =axis ;
    %
    % AAA(2)  =NMODES_SHOW ;
    % axis(AAA) ;
    
    figure(nfigure)
    hold on
    
    %dbstop('32')
    h = plot([ 1:length(S)],[ normEf2(2:end); eSVD]/nTOTAL*100,[COLOR,'-*']');
    xlabel('MODES')
    ylabel([' SVD error (%)' ])
    if ~isempty(NMODES_SHOW)
        AAA =axis ;
        AAA(2)  =NMODES_SHOW ;
        axis(AAA) ;
    end
    
    
    figure(nfigure+1000)
    hold on
    h2 = plot([ 1:length(S)],[ log10([normEf2(2:end); eSVD]/nTOTAL)],[COLOR,'-*']');
    xlabel('MODES')
    ylabel([' logarithm SVD error' ])
    if ~isempty(NMODES_SHOW)
        AAA =axis ;
        AAA(2)  =NMODES_SHOW ;
        axis(AAA) ;
    end
end
