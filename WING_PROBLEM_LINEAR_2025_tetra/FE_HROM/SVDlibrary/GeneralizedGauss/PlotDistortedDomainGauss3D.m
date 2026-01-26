function PlotDistortedDomainGauss3D(DATAOUT_irreg_gauss,DATAOUT_gauss,wGAUSS,xGAUSS,VSinv,DATALOC,xMAT,errorGAUSS,w_gauss_q)

if nargin == 0
    load('tmp1.mat')
end

if  ~isempty(DATAOUT_irreg_gauss)
    % Plot points in deformed space
    %
    figure(89)
    hold on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Distorted domain ')
    
    % Plotting the boundary of the distorted domain 
    % --------------------------
    xLIM = DATALOC.xLIM ; 
    nPLOC = 100 ; 
    iDIM = 1; 
    xxVAR = linspace(xLIM(iDIM,1),xLIM(iDIM,2),nPLOC)' ; 
    iDIM = 2; 
    yyVAR = linspace(xLIM(iDIM,1),xLIM(iDIM,2),nPLOC)' ;
    iDIM = 3; 
    zzVAR = linspace(xLIM(iDIM,1),xLIM(iDIM,2),nPLOC)' ;
   
    xxFIXmin =  xLIM(1,1)*ones(size(xxVAR)) ; 
    xxFIXmax =  xLIM(1,2)*ones(size(xxVAR)) ; 
    
    yyFIXmin =  xLIM(2,1)*ones(size(yyVAR)) ; 
    yyFIXmax =  xLIM(2,2)*ones(size(yyVAR)) ; 
    
     zzFIXmin =  xLIM(3,1)*ones(size(zzVAR)) ; 
    zzFIXmax =  xLIM(3,2)*ones(size(zzVAR)) ; 
    
 
    
    xBORDER = [xxVAR,yyFIXmin,zzFIXmin
              xxFIXmax,yyVAR,zzFIXmin
              xxVAR(end:-1:1),yyFIXmax,zzFIXmin
              xxFIXmin,yyVAR(end:-1:1),zzFIXmin 
              %
              xxVAR,yyFIXmin,zzFIXmax
              xxFIXmax,yyVAR,zzFIXmax
              xxVAR(end:-1:1),yyFIXmax,zzFIXmax
              xxFIXmin,yyVAR(end:-1:1),zzFIXmax
              %
              xxFIXmin,yyFIXmin,zzVAR
              xxFIXmax,yyFIXmin,zzVAR
              xxFIXmax,yyFIXmax,zzVAR 
               xxFIXmin,yyFIXmax,zzVAR ] ; 
          
          
          
          
    DATALOC.EVALUATE_SHAPE_FUNCTIONS = 1; 
    EVALUATE_GRADIENT = 0 ; 
    EVALUATE_FUN = 0 ; 
          
    [PHIk_y dPHIk_y,DATAOUTloc]= ...
    EvaluateBasisFunctionAtX(xBORDER,VSinv,DATALOC,EVALUATE_GRADIENT,EVALUATE_FUN) ; 
    
    
     plot3(DATAOUTloc.COORdistorted(:,1),DATAOUTloc.COORdistorted(:,2),DATAOUTloc.COORdistorted(:,3),'r.')  ;
    
    
    
    
    % Shape functions (our points)
    
    
    
       [PHIk_y dPHIk_y,DATAOUTloc_mine]= ...
    EvaluateBasisFunctionAtX(xGAUSS,VSinv,DATALOC,EVALUATE_GRADIENT,EVALUATE_FUN) ; 
    
    
    Xtransf = DATAOUTloc_mine.COORdistorted(:,1) ;
    Ytransf =  DATAOUTloc_mine.COORdistorted(:,2) ;
    Ztransf =  DATAOUTloc_mine.COORdistorted(:,3) ;
    
    hOUR = plot3(Xtransf,Ytransf,Ztransf,'go') ;
    
         [PHIk_y dPHIk_y,DATAOUTloc_gauss]= ...
    EvaluateBasisFunctionAtX(xMAT,VSinv,DATALOC,EVALUATE_GRADIENT,EVALUATE_FUN) ; 
    
    
    Xtransf_gauss =  DATAOUTloc_gauss.COORdistorted(:,1) ;
    Ytransf_gauss =  DATAOUTloc_gauss.COORdistorted(:,2) ;
    Ztransf_gauss =  DATAOUTloc_gauss.COORdistorted(:,3) ;
    
    h = plot3(Xtransf_gauss,Ytransf_gauss,Ztransf_gauss,'bo') ;
    legend(h,{})
    
    
    legend([hOUR,h],{['Our rule (',num2str(length(Xtransf)),' points)'],...
        ['Gauss rule (',num2str(length(Xtransf_gauss)),' points, error = ',num2str(errorGAUSS),' %)']})
    
    SHOW_TEXT = 0; 
    if SHOW_TEXT ==1
    for ipoints = 1:length(wGAUSS)
        w = wGAUSS(ipoints)*DATAOUTloc_mine.detJ(ipoints) ;
        text(Xtransf(ipoints),Ytransf(ipoints),num2str(w)) ;
    end
    
    
    for ipoints = 1:length(w_gauss_q)
        w = w_gauss_q(ipoints)*DATAOUTloc_gauss.detJ(ipoints) ;
        text(Xtransf_gauss(ipoints),Ytransf_gauss(ipoints),num2str(w)) ;
    end
    end
end




