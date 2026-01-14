function PlotStressPoints(iply,Fenv,stresses)

   ifigure = 100+iply ; 
      figure(ifigure)
        hold on
        title(['STRESSES -->PLY ',num2str(iply)])
        %
        % Points ellipsoid 
        LEGENDS  = {'\sigma_1','\sigma_2','\tau_{12}'} ;
         sELIP = [1 0 0 1 1 0  
                  0 1 0 1 0 1
                  0 0 1 0 1 1]; 
          sELIP = [sELIP -sELIP];
                      
         F_stress = Fenv*sELIP ;
         FailureCrit = sqrt(sum(sELIP.*F_stress,1)) ;
         
         
         sELIP = bsxfun(@times,sELIP',1./FailureCrit')' ; 
        %
        DATAPLOT.PLOT_POINTS_DEFINING_SURF = 0 ; 
        [Fenv radii_stress evecs_stress chi2_stress] = PlotEnvelopeFailure(sELIP(1,:)',sELIP(2,:)',sELIP(3,:)',...
            LEGENDS,ifigure) ; 
        
        % Now we plot all points 
        % ..................... .
        
        hold on
        plot3(stresses(1,:),stresses(2,:),stresses(3,:),'*')
     