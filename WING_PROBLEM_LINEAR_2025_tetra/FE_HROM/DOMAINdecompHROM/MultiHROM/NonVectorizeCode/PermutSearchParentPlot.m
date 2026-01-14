function PermutSearchParentPlot(XeRELs,XrefAUG,PLOTTR,XeRELsROTATED)

if nargin == 0
    load('tmp1.mat')
    PLOTTR = 1; 
%     ANGLE =- 180 ; 
%     Q = [cos(ANGLE) -sin(ANGLE)  0
%          sin(ANGLE) cos(ANGLE)  0
%          0  0  1] ; 
%      
%     XeRELsROTATED = Q'*XeRELs ; 
end

if PLOTTR == 1
    
    
    if size(XeRELs,1) == 2
        
        figure(124)
        hold on
        xlabel('x')
        ylabel('y')
        XeAUG = [XeRELs,XeRELs(:,1)] ;
        
        h1 = plot(XrefAUG(1,:),XrefAUG(2,:)) ;
        h2= plot(XeAUG(1,:),XeAUG(2,:)) ;
        for ipoints = 1:size(XeRELs,2)
            text(XrefAUG(1,ipoints),XrefAUG(2,ipoints),['R',num2str(ipoints)])
        end
        for ipoints = 1:size(XeRELs,2)
            text(XeRELs(1,ipoints),XeRELs(2,ipoints),['     E',num2str(ipoints)]);
        end
        
    else
        
        
        
        % HEXAHEDRA ELEMENTS
        
        SHOW_ROTATED= 0;
        if SHOW_ROTATED == 1
            XeRELs =  XeRELsROTATED ;
        end
        
        figure(125)
        hold on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        lines_PLOT = [1,2; 2 3; 3 4; 4 1;  5 6; 6 7; 7 8; 8 5; 1 5;  2  6; 3 7 ; 4 8  ]; 
        %XeAUG = [XeRELs,XeRELs(:,1)] ;
        
        for iplot = 1:length(lines_PLOT)
            ind1 = lines_PLOT(iplot,1) ; ind2 = lines_PLOT(iplot,2) ;
            xx =  [XrefAUG(1,ind1),XrefAUG(1,ind2) ] ;
            yy = [XrefAUG(2,ind1),XrefAUG(2,ind2) ] ;
            zz = [XrefAUG(3,ind1),XrefAUG(3,ind2) ] ;
            plot3(xx,yy,zz) ; 
        end
        for ipoints = 1:8
            text(XrefAUG(1,ipoints),XrefAUG(2,ipoints),XrefAUG(3,ipoints),['  R',num2str(ipoints)])
        end
        
        
        for iplot = 1:length(lines_PLOT)
            ind1 = lines_PLOT(iplot,1) ; ind2 = lines_PLOT(iplot,2) ;
            xx =  [XeRELs(1,ind1),XeRELs(1,ind2) ] ;
            yy = [XeRELs(2,ind1),XeRELs(2,ind2) ] ;
            zz = [XeRELs(3,ind1),XeRELs(3,ind2) ] ;
            plot3(xx,yy,zz);
        end
        for ipoints = 1:8
            text(XeRELs(1,ipoints),XeRELs(2,ipoints),XeRELs(3,ipoints),['   (X',num2str(ipoints),')']);
        end
        
        
        
    end
    
end