function  [U,S,V,h,h2,Ui_Si] = SVD_dom(SNAP,nfigure,LEGENDG,NMODES_SHOW,COLOR,DATA,...
    INDICES_PROJECTS)


if nargin ==  0
    load('tmp1.mat')
end

h= [] ; h2 = [] ;
DATA = DefaultField(DATA,'SVD_RANDOM',0) ;
DATA = DefaultField(DATA,'TOL_LOC',0) ;
DATA = DefaultField(DATA,'TOLERANCE_SVD_PROJECT',[] ) ;
DATA = DefaultField(DATA,'TOLERANCE_SVD_SLICE',cell(size(SNAP)) ) ;

NMODES = DATA.NMODES_PROJECT_LOC; %(INDICES_PROJECTS) ;
if iscell(NMODES)
    for iproj = 1:length(NMODES)
        nmodesMAX = size(SNAP{iproj},2) ;
        if ~isempty(NMODES{iproj})
            % NMODES{iproj} = nmodesMAX ;
            %else
            %    NMODES{iproj} = min(nmodesMAX,NMODES{iproj})  ;
        end
    end
end
%NMODES = cell2mat(NMODES) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
DATA = DefaultField(DATA,'PARTITION_COLUMNS',[])  ; 
if ~isempty(DATA.PARTITION_COLUMNS)
    NSTEPS = DATA.PARTITION_COLUMNS(INDICES_PROJECTS) ;  % For dealing with multiple time steps
    TOLERANCE_SVD_SLICE = DATA.TOLERANCE_SVD_SLICE(INDICES_PROJECTS) ;
    % Proceed to split SNAP into submatrices
    % --------------------------------------
    for iproj = 1:length(SNAP)
        if NSTEPS(iproj) >2
          
            
           
            if ~isempty(TOLERANCE_SVD_SLICE{iproj})
                [nnn,mmm] = size(SNAP{iproj});
                 SNAP{iproj} = mat2cell(SNAP{iproj},nnn,COLS) ;
                %  Number of domains
                ndom = mmm/NSTEPS(iproj) ;
                COLS = NSTEPS(iproj)*ones(1,ndom) ;
                DATA.TOL_BLOCK_SVD{iproj} = TOLERANCE_SVD_SLICE{iproj}*ones(size(SNAP{iproj})) ;
                
            end
        end
    end
    
    
else
    % If TOLERANCE_SVD_SLICE{iproj}
    
   DATA.TOL_BLOCK_SVD= DATA.TOLERANCE_SVD_SLICE(INDICES_PROJECTS)   ;
   for iblock = 1:length(DATA.TOL_BLOCK_SVD)
       if isempty(DATA.TOL_BLOCK_SVD{iblock})
           DATA.TOL_BLOCK_SVD{iblock} = 0 ; 
       end
   end
 % Tolerance for each block(project) 
 
end

DATALOC.TOL_BLOCK_SVD = DATA.TOL_BLOCK_SVD ;
DATALOC.RELATIVE_SVD = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA = DefaultField(DATA,'TOL_SVD_GLO',0) ; 
DATALOC.TOL_SVD_GLO = DATA.TOL_SVD_GLO ; 

[U,S,V,Ui_Si] = RSVDcol(SNAP,NMODES,DATALOC) ;
eSVD = 0 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% elseif ~isempty(DATA.TOLERANCE_SVD_PROJECT)
%
%     [U,S,V,eSVD] = SVD_accuracy_project(SNAP,DATA,COLUMNS_RVEloc) ;
%
%
% else
%     error('Option not implemented')
%
% end

DATA = DefaultField(DATA,'PLOT_RIGHT_SINGULAR_VECTORS',0); % .PLOT_RIGHT_SINGULAR_VECTORS =1 ;

if DATA.PLOT_RIGHT_SINGULAR_VECTORS >0 &&  ~isempty(S)
    NMODES_RIGHT = min(DATA.PLOT_RIGHT_SINGULAR_VECTORS,length(S)) ;
    LLL = DATA.LEGEND_GRAPHS ;
    LLLeg = {} ;
    %     if DATA.LEGEND_GRAPHS(1) == 'D'
    %         nfigures = 400 ;
    %     else
    %         nfigures  = 401 ;
    %     end
    figure(nfigure.RIGHTSV)
    hold on
    xlabel('Domain')
    ylabel('Right S. Val.')
    colores = ColoresMatrix(NMODES_RIGHT );
    title(['Right singular v. x Sing. Val., ',DATA.LEGEND_GRAPHS])
    
    Snorm = sqrt(sum(S.^2)) ; 
    
    for iii = 1:NMODES_RIGHT
        hLOC(iii) = plot(V(:,iii)*S(iii)/Snorm,'Color',colores(iii,:)) ;
        LLLeg{iii} = ['i=',num2str(iii),',  s =',num2str(S(iii)/Snorm,4)] ;
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
    
    figure(nfigure.ERROR)
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
    
    
    figure(nfigure.ERROR+1)
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
