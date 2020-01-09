function SuperEllipsePlaneRhoTxi
mxMax = 0.99;
myMax = 0.99;
mxMin = 0;
myMin = 0;
c = SuperEllipseParamsRelator.c);
qV = [2 10^6];
color = {'b','r'};
outFile = '/home/alex/Dropbox/PaperStress/AddmissibleSpace/';

% for iq = 1:length(qV)
%     q = qV(iq);
%     rhoMin = 1 - c(q)*mxMax*myMax;
%     rho(:,1) = linspace(rhoMin,1, 100);
%     txiUB(:,1) = atan(mxMax^2*c(q)./(1-rho));
%     txiLB(:,1) = atan((1-rho)/(myMax^2*c(q)));
%     plot(rho,[txiLB txiUB],['-',color{iq}])
%     set(gca,'ytick',[0:pi/8:pi/2]) % where to set the tick marks
%     set(gca,'yticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
%     xlabel('\rho')
% end

myCase = 'A';

switch myCase
    case 'A'
        
        %figure(1)
        for iq = 1:length(qV)
            q = qV(iq);
            txi = linspace(0,pi/2,1000);
            s.q = q;
            s.txi = txi;
            s.mxMin = mxMin;
            s.myMin = myMin;  
            s.mxMax = mxMax;            
            s.myMax = myMax;  
            rhoBounds = SuperEllipseRhoBoundsComputer(s);
            [rhoMin,rhoMax] = rhoBounds.compute();   
            
            figure()
            h1 = plot(txi,rhoMin,['-',color{iq}]);
            set(h1,'LineWidth',2);        
            hold on
            h2 = plot(txi,rhoMax,['-',color{iq}]);
            set(h2,'LineWidth',2);
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})            
            xlabel('\xi')
            ylabel('\rho')
            set(gca,'FontSize',18) 
            axis([0 pi/2 0 1])
            print([outFile,'AdmissibleSpaceForQ',num2str(q)],'-dpdf')
        end
    case 'B'
        figure(1)
        hold on
        for iq = 1:length(qV)
            q = qV(iq);
            txi = linspace(0,pi/2,1000);
            rhoMinMx = 1- mxMax^2*c(q)./tan(txi);
            rhoMinMy = 1- myMax^2*c(q)*tan(txi);
            rhoMin = max(rhoMinMx,rhoMinMy);
            h1 = plot(txi,rhoMin,['-',color{iq}]);
            set(h1,'LineWidth',2);        
            set(gca,'xtick',[0:pi/8:pi/2]) % where to set the tick marks
            set(gca,'xticklabels',{'0','\pi/8','\pi/4','3\pi/8','\pi/2'})
            xlabel('\xi')
            ylabel('\rho')
            set(gca,'FontSize',18)             
            axis([0 pi/2 0 1])
            print([outFile,'AdmissibleSpace',num2str(q)],'-dpdf')
        end
end


end