classdef ProposedOptimalSuperEllipse < handle

    properties (Access = private)
       mxV
       myV
       qV
       alpha
       beta
       outputPath
    end
    
    methods (Access = public)
        
        function obj = ProposedOptimalSuperEllipse()
            obj.outputPath = '/home/alex/Dropbox/PaperStress/';
            obj.obtainOptimalSuperEllipseValues();
            obj.obtainOptimalCoeficientsForAnalyticalSuperEllipse(); 
            obj.plotAnalyticalSuperEllipse();
            obj.plotAnalyticalQvsMaxMxMy();
        end        
        
        function obtainOptimalSuperEllipseValues(obj)
            superEllipse = PonderatedOptimalSuperEllipseComputer();
            superEllipse.compute();
            obj.mxV = superEllipse.mxV;
            obj.myV = superEllipse.myV;
            obj.qV  = superEllipse.qMean;            
        end
        
        function obtainOptimalCoeficientsForAnalyticalSuperEllipse(obj)
            qAnalytical = @(alpha,beta) obj.analyticalExponent(alpha,beta,obj.mxV,obj.myV);
            options = optimoptions('fmincon','Display','iter');
            problem.options = options;
            problem.solver = 'fmincon';
            problem.objective = @(x) norm(qAnalytical(x(1),x(2)) - obj.qV)/norm(obj.qV);
            problem.x0 = [6,20];
            problem.ub = [50,50];
            problem.lb = [0,0];
            x = fmincon(problem);
            obj.alpha = x(1);
            obj.beta = x(2);
        end
        
        function plotAnalyticalSuperEllipse(obj)
            f = obj.plotExponent(obj.alpha,obj.beta);
            p = surfPrinter(f);
            xlabel('$m_1$','Interpreter','latex');
            ylabel('$m_2$','Interpreter','latex');
            zlabel('$q$','Interpreter','latex');
            fName = fullfile(obj.outputPath,'AnalyticalOptimalSuperEllipse');
            p.print(fName)            
        end
        
        function plotAnalyticalQvsMaxMxMy(obj)
            mx = linspace(0.01,0.99,100);
            my = linspace(0.01,0.99,100);
            q = obj.analyticalExponent(obj.alpha,obj.beta,mx,my);
            f = figure();
            h{1} = plot(mx,q,'-+');
            xlabel('$\max(m_1,m_2)$','Interpreter','latex');
            ylabel('$q$','Interpreter','latex');
            fName = fullfile(obj.outputPath,'AnalyticalQvsMaxMxMy');
            p = plotPrinter(f,h);
            p.print(fName);            
        end
        
        function f = plotExponent(obj,alpha,beta)
            f = figure();
            nMx = 100;
            nMy = 100;           
            mxv = linspace(0.01,0.99,nMx);
            myv = linspace(0.01,0.99,nMy);
            q = zeros(nMx,nMy);
            for ix = 1:nMx
                for iy = 1:nMy
                    mx = mxv(ix);
                    my = myv(iy);
                    q(ix,iy) = obj.analyticalExponent(alpha,beta,mx,my);
                end
            end
            surf(mxv,myv,q')           
        end
        
    end
    
    methods (Access = private, Static)
        
        function q = analyticalExponent(alpha,beta,m1,m2)
            qmin = 2;
            qmax = 32;
            q = max(qmin,min((1./(1-max(m1,m2).^beta)).^alpha,qmax));
        end
        
    end
    
end