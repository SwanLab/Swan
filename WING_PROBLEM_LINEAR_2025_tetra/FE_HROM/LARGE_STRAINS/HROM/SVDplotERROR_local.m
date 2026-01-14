function SVDplotERROR_local(S,ifig)


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
    %    ifig = 3000;
        figure(ifig)
        hold on
        title(['SVD error (%), internal work snapshots ']) ;
        eSVD = 0 ;
        %dbstop('32')
        
        
        normREL = normEf2/nTOTAL;
     %   disp(['ERROR SVD (over 1) icluster = ',num2str(icluster),' var =',nameVAR])
     %   normREL
        
        h = plot([ 1:length(S)],[ normEf2(2:end); eSVD]/nTOTAL*100 );
        xlabel('MODES')
        ylabel([' SVD error (%)' ])