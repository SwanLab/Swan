classdef ComparingSuperEllipseVsRectangle < handle
       
    properties (Access = private)
       mesh 
       dataRes
       m1S
       m2S
       m1R
       m2R
    end
    
    properties (Access = private)
        path
        fileNameMesh
        fileResName
        pNorm
        microCase
    end
    
    methods (Access = public)
        
        function obj = ComparingSuperEllipseVsRectangle()
            obj.init()
            obj.wrapResAndMesh();
            obj.computeDesignParams();
            obj.computeStresses();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.microCase = 'SuperEllipse';'Rectangle';
            meshType  = 'Medium';'Small';
            finalIter = '589';'485';'405';'250';'296';
            fCase = [obj.microCase,'Rotation',meshType];
            obj.path = ['/media/alex/My Passport/LatticeResults/StressNorm',fCase,'/'];
            obj.fileResName  = ['ExperimentingPlot',finalIter];
            obj.fileNameMesh = ['CantileverSquare',meshType];
            obj.pNorm = 2:2:32;
        end
        
        function wrapResAndMesh(obj)
            s.folderPath = obj.path;
            s.fileName   = obj.fileResName;
            w = WrapperMshResFiles(s);
            w.compute();
            obj.mesh    = w.mesh;
            obj.dataRes = w.dataRes;
        end
        
        function computeDesignParams(obj)
            switch obj.microCase
                case 'Rectangle'
                   obj.computeDesignParamsFromrectangleExample();
                case 'SuperEllipse'
                   obj.computeDesignParamsFromSuperEllipseExample();
            end
        end
        
        function computeDesignParamsFromSuperEllipseExample(obj)
            obj.m1S = obj.dataRes.DesignVar1;
            obj.m2S = obj.dataRes.DesignVar2;
            [obj.m1R,obj.m2R] = obj.computeRectangleM1M2FromSuperEllipseM1M2(obj.m1S,obj.m2S);
        end
        
        function computeDesignParamsFromrectangleExample(obj)
            obj.m1S = obj.dataRes.DesignVar1;
            obj.m2S = obj.dataRes.DesignVar2;
            [obj.m1R,obj.m2R] = obj.computeRectangleM1M2FromSuperEllipseM1M2(obj.m1S,obj.m2S);
        end
                
        function computeStresses(obj)
           valueR = obj.computeStressNormRectangle();
           valueS = obj.computeStressNormSuperEllipse();
           f = figure();
           hold on
           p{1} = plot(obj.pNorm,valueR);
           p{2} = plot(obj.pNorm,valueS);
           legend({'Rectangle','SuperEllipse'},'Location','Best');
           xlabel('p');
           ylabel('||\sigma||_p');
           pP = plotPrinter(f,p);
           fName = '/home/alex/Dropbox/GregMeeting30Octubre/Comparison';
           pP.print(fName)
        end
        
        function value = computeStressNormSuperEllipse(obj)
            s.m1            = obj.m1S;
            s.m2            = obj.m2S;
            s.alpha         = obj.dataRes.AlphaGauss';
            s.vademecumName = 'SuperEllipseQOptAnalytic';
            s.mesh          = obj.mesh;
            s.fileName      = obj.fileNameMesh;
            s.pNorm         = obj.pNorm;
            sC = StressNormFromVademecumComputer(s);
            value = sC.compute();
        end        
        
        function value = computeStressNormRectangle(obj)
            s.m1            = obj.m1R;
            s.m2            = obj.m2R;
            s.alpha         = obj.dataRes.AlphaGauss';
            s.vademecumName = 'SuperEllipseQMax';
            s.mesh          = obj.mesh;
            s.fileName      = obj.fileNameMesh;
            s.pNorm         = obj.pNorm;
            sC = StressNormFromVademecumComputer(s);
            value = sC.compute();
        end
        
    end
    
    methods (Access = private, Static)
       
        function [m1R,m2R] = computeRectangleM1M2FromSuperEllipseM1M2(m1S,m2S)           
            sE   = SuperEllipseParamsRelator;
            s.type  = 'Optimal';
            s.m1    = m1S;
            s.m2    = m2S;
            sM = SmoothingExponentComputer.create(s);
            qS(:,1) = sM.compute();
            rho  = sE.rho(m1S,m2S,qS);
            xi   = sE.xi(m1S,m2S);
            s.type  = 'Given';
            s.q    = 32*ones(size(m1S));
            sM = SmoothingExponentComputer.create(s);
            qR(:,1) = sM.compute();
            m1R = sE.mx(xi,rho,qR);
            m2R = sE.my(xi,rho,qR);
        end
        
        function [m1S,m2S] = computeSuperEllipseM1M2FromRectangleM1M2(m1R,m2R)
            sE   = SuperEllipseParamsRelator;
            rho  = sE.rho(m1R,m2R,q);
            xi   = sE.xi(m1R,m2R);
            s.type  = 'Optimal';
            s.m1    = [];
            s.m2    = [];
            sM = SmoothingExponentComputer.create(s);
            qS(:,1) = sM.computeQ(rho,xi);
            m1S = sE.mx(xi,rho,qS);
            m2S = sE.my(xi,rho,qS);
        end                
        
    end
    
end