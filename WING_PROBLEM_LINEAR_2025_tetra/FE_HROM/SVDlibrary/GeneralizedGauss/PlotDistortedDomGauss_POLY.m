function PlotDistortedDomGauss_POLY(DATAOUT_irreg_gauss,DATAOUT_gauss,xGAUSS,wGAUSS,w_gauss_q,xMAT,errorGAUSS)




if  ~isempty(DATAOUT_irreg_gauss)
    % Plot points in deformed space
    %
    figure(89)
    hold on
    xlabel('x')
    ylabel('y')
    title('Distorted domain ')
    
    X = [DATAOUT_irreg_gauss.Xcoor,DATAOUT_irreg_gauss.Xcoor(:,1)] ;
    plot(X(1,:),X(2,:),'r')  ;
    % Shape functions (our points)
    
    N = DATAOUT_gauss.Nshape(xGAUSS(:,1),xGAUSS(:,2))  ;
    
    Xtransf = N*DATAOUT_gauss.Xcoor(1,:)' ;
    Ytransf = N*DATAOUT_gauss.Xcoor(2,:)' ;
    
    hOUR = plot(Xtransf,Ytransf,'rx') ;
    
    N = DATAOUT_gauss.Nshape(xMAT(:,1),xMAT(:,2))  ;
    
    
    Xtransf_gauss = N*DATAOUT_irreg_gauss.Xcoor(1,:)' ;
    Ytransf_gauss = N*DATAOUT_irreg_gauss.Xcoor(2,:)' ;
    
    h = plot(Xtransf_gauss,Ytransf_gauss,'bo') ;
    legend(h,{})
    
    
    legend([hOUR,h],{['Our rule (',num2str(length(Xtransf)),' points)'],...
        ['Gauss rule (',num2str(length(Xtransf_gauss)),' points, error = ',num2str(errorGAUSS),' %)']})
    
    
    for ipoints = 1:length(wGAUSS)
        w = wGAUSS(ipoints)*DATAOUT_gauss.detJ(ipoints) ;
        text(Xtransf(ipoints),Ytransf(ipoints),num2str(w)) ;
    end
    
    
    for ipoints = 1:length(w_gauss_q)
        w = w_gauss_q(ipoints)*DATAOUT_gauss.detJ(ipoints) ;
        text(Xtransf_gauss(ipoints),Ytransf_gauss(ipoints),num2str(w)) ;
    end
end
