function PlotSVDmodesFitting(DATALOC,iMODE,S,Ux,Vy,Ux_splineCOEFF,Vy_splineCOEFF,x,y)

PlotVector = 0 ; 
xmin = min(x); 
xmax = max(x) ; 

ymin = min(y) ; 
ymax = max(y) ;
mPOINTSx = 1000 ; 
mPOINTSy = 500 ; 
xFINE = linspace(xmin,xmax,mPOINTSx) ; 

yFINE = linspace(ymin,ymax,mPOINTSy) ; 

     figure(100+iMODE)
     % x functions 
     % -------------
     subplot(2,1,1)
     %
     hold on
     h = [] ;
     xlabel('x')
     ylabel('U(x)')
     title(['Left singular vectors, ', 'FUN =',num2str(iMODE)])
     LLL = {} ;
     for i = 1:length(S)
         
         if PlotVector ==1 
         h(end+1) = plot(x,Ux(:,i)) ;
         LLL{end+1} = ['U_',num2str(i),'; s = ',num2str(S(i))] ;
         end
         %         
         % Approximated
       %  if PlotApproximated == 1
         Ux_approx = ppval(Ux_splineCOEFF{i},xFINE) ;
         h(end+1) = plot(xFINE,Ux_approx) ;
         LLL{end+1} = ['U_',num2str(i),'; s = ',num2str(S(i)),'(spl)'  ] ;     
       %  end
         
     end  
     
         legend(h,LLL)
     
     % y functions 
     % -------------
    
     % -------------
     subplot(2,1,2)
     %
     hold on
     h = [] ;
     xlabel('y')
     ylabel('V(x)')
     title(['Right singular vectors, ', 'FUN =',num2str(iMODE)])
     LLL = {} ;
     for i = 1:length(S)
          if PlotVector ==1 
         h(end+1) = plot(y,Vy(:,i)) ;
         LLL{end+1} = ['V_',num2str(i),'; s = ',num2str(S(i))] ;
          end
         %         
         % Approximated
     %    if PlotApproximated == 1 
         Vy_approx = ppval(Vy_splineCOEFF{i},yFINE) ;
         h(end+1) = plot(yFINE,Vy_approx) ;
         LLL{end+1} = ['V_',num2str(i),'; s = ',num2str(S(i)),' (spl)'  ] ;         
      %   end
     end  
     
         legend(h,LLL)
 