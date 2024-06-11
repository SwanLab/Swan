classdef RegularizedSuperEllipseExponentPlotter < handle
    
    properties (Access = private)
        vademecum
        mesh
        errorF
        qAveraged
        qRegularized
        qComputer        
    end
    
    methods (Access = public)
        
        function obj = RegularizedSuperEllipseExponentPlotter()
            obj.createQcomputer();
            obj.loadVademecum();
            obj.createMesh();
            obj.createErrorFunction();
            obj.computeNumericalAveragedValue();
            obj.computeQregularized();
            obj.plotQforSymmetricHoles();
            obj.plotSuperEllipseOptimalExponent();
        end
        
    end
    
    methods (Access = private)
                
        function createQcomputer(obj)
            s.m1 = [];
            s.m2 = [];
            obj.qComputer = SmoothingExponentComputerOptimal(s);            
        end        
                
        function loadVademecum(obj)
            obj.vademecum = VademecumReader();            
        end
        
        function createMesh(obj)
            xi = obj.vademecum.xiV;
            rho = obj.vademecum.rhoV;
            obj.mesh = obj.obtainMesh(xi,rho);            
        end
        
        function computeNumericalAveragedValue(obj)
            s.vademecum = obj.vademecum;
            p = PonderatedOptimalSuperEllipseComputer(s);
            p.compute();
            obj.qAveraged = p.qMean;            
        end
        
        function error = computeError(obj,a,b,r,qS)
            q = obj.computeQ(a,b,r,qS); 
            error = obj.errorF(q);
        end
        
        function q = computeQ(obj,a,b,r,q)
             obj.qComputer.setParamsValues(a,b,r,q)  
             xi  = obj.vademecum.xiV;
             rho = obj.vademecum.rhoV;             
             q = obj.qComputer.computeQ(xi,rho);                
        end
        
        function [alpha,beta,rhoQmin,qSoptMin] = computeBestParameters(obj)
            s.errorComputer = @(a,b,c,d) obj.computeError(a,b,c,d);
            r = RegularizedExponentBestCoeffComputer(s);
            [alpha,beta,rhoQmin,qSoptMin] = r.compute();
        end
                
        function createErrorFunction(obj)
            s.type              = 'SIMPLE';
            s.mesh              = obj.mesh;
            s.backgroundMesh    = obj.mesh;
            s.globalConnec      = obj.mesh.connec;
            s.npnod = obj.mesh.nnodes;
            int = Integrator.create(s);
            int.computeLHS();
            p = 2;            
            normF = @(x) ((int.computeL2norm(x)).^(p))^(1/p);            
            obj.errorF = @(x) normF(abs(x - obj.qAveraged))/normF(obj.qAveraged);
        end        
        
        function computeQregularized(obj)
            [alpha,beta,rhoQmin,qSoptMin] = obj.computeBestParameters();                       
            q = obj.computeQ(alpha,beta,rhoQmin,qSoptMin);            
            obj.qRegularized = q;
        end
        
        function plotQforSymmetricHoles(obj)
            xi     = obj.vademecum.xiV;
            isSym  = abs(xi-pi/4) < 0.01;
            isSym  = isSym(2:end);
            rhoSym = obj.vademecum.rhoV(isSym);
            qAsym  = obj.qAveraged(isSym);
            qRsym  = obj.qRegularized(isSym);
            f = figure();
            hold on
            h{1} = plot(rhoSym,qAsym,'-+');
            h{2} = plot(rhoSym,qRsym,'-+');
            xlabel('$\rho$','Interpreter','Latex')
            ylabel('$q(\xi = \pi/4,\rho)$','Interpreter','Latex')            
            legN = '$\textrm{Numerical smoothing exponent} \, q_N$';
            legA = '$\textrm{Analytical smoothing exponent} \, q_A$';            
            legend({legN,legA},'Interpreter','Latex','Location','Best');            
            outPutPath = '/home/alex/git-repos/MicroStructurePaper/';
            outputName = [outPutPath,'qMaxM1M2AnalyticalNumerical'];
            printer = plotPrinter(f,h);
            printer.print(outputName);            
        end
        
        function plotSuperEllipseOptimalExponent(obj)
            s.title = 'Proposed ';
            s.fileName = '/home/alex/git-repos/MicroStructurePaper/ProposedOptimal';
            s.rhoV = rho;
            s.xiV = xi;
            s.value = qA;
            s.qMean = qA;
            p = SuperEllipseExponentPlotter(s);
            p.plot();            
        end
        
        function m = obtainMesh(obj,x,y)
            s.coord = [x,y];
            s.connec = obj.obtainConnec(x,y);
            m = Mesh.create(s);
        end
        
        function connec = obtainConnec(obj,x,y)
            connec = delaunay(x,y);
            connec = obj.obtainQualityElements(connec,x,y);
        end        
        
    end
    
    methods (Access = private, Static)
               
        function connec = obtainQualityElements(connec,x,y)
            s.coord = [x,y];
            s.connec = connec;
            m = Mesh.create(s);
            qua = m.computeElementQuality';
            isQ = qua > 0.02;
            connec = connec(isQ,:);
        end
        
    end
end