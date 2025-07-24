classdef VademecumCellVariablesCalculator < handle
    
    properties (Access = public)
       variables 
    end
    
    properties (Access = private)
        fileNames        
        homog        
        iMxIndex
        iMyIndex
        mxV
        myV
        iter
        freeFemSettings
        print
        smoothingExponentSettings
    end
    
    methods (Access = public)
        
        function obj = VademecumCellVariablesCalculator(d)
            obj.init(d);
        end
        
        function computeVademecumData(obj)
            nMx = length(obj.mxV);
            nMy = length(obj.myV);
            for imx = 1:nMx
                for imy = 1:nMy
                    obj.storeIndex(imx,imy);
                    obj.iter = (imy + nMx*(imx -1));
                    disp([num2str(obj.iter/(nMx*nMy)*100),'% done']);                    
                   % obj.generateMeshFile();
                    obj.computeNumericalHomogenizer();
                    obj.obtainHomogenizerData();
                end
            end
        end        
        
        function saveVademecumData(obj)
           d = obj.getData();
           fN = obj.fileNames.fileName;
           pD = obj.fileNames.outPutPath;
           file2SaveName = [pD,'/',fN,'.mat'];
           save(file2SaveName,'d');
        end
        
        function d = getData(obj)
            d.variables        = obj.variables;
            d.domVariables.mxV = obj.mxV;
            d.domVariables.myV = obj.myV;
            d.outPutPath    = obj.fileNames.outPutPath;
            d.fileName      = obj.fileNames.fileName;
        end
         
    end    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.computeFileNames(d);
            nMx = d.nMx;
            nMy = d.nMy;
            obj.mxV = linspace(d.mxMin,d.mxMax,nMx);
            obj.myV = linspace(d.myMin,d.myMax,nMy);
            obj.print = d.print;
            obj.smoothingExponentSettings = d.smoothingExponentSettings;
        end
        
        function computeFileNames(obj,d)
            fN = d.fileName;
            oP = d.outPutPath;
            pD = fullfile(pwd,'Output',fN);
            fNs.fileName    = fN;
            fNs.outPutPath  = oP;
            fNs.printingDir = pD;                        
            obj.fileNames = fNs;
        end        
        
        function storeIndex(obj,imx,imy)
            obj.iMxIndex = imx;
            obj.iMyIndex = imy;
        end
        
        
        function q = computeCornerSmoothingExponent(obj)
            s = obj.smoothingExponentSettings;
            s.m1 = obj.mxV(obj.iMxIndex);
            s.m2 = obj.myV(obj.iMyIndex);
            qComputer = SmoothingExponentComputer.create(s);
            q = qComputer.compute();
        end
             
        function computeNumericalHomogenizer(obj)
%             d.gmsFile = obj.fileNames.gmsFile;
%             d.outFile = [obj.fileNames.fileName];
%             d.print   = obj.print;
%             d.iter = obj.iter;
%             nH = NumericalHomogenizerCreatorFromGmsFile(d);
%             obj.homog = nH.getHomogenizer();
            
            s.fileName = ['OptimaSuperEllipseCase',num2str(obj.iter)];
            
            mx = obj.mxV(obj.iMxIndex);
            my = obj.mxV(obj.iMyIndex);           
            q = obj.computeCornerSmoothingExponent();
            
            %s.rho   = SuperEllipseParamsRelator.rho(mx,my,q);
            %s.txi   = SuperEllipseParamsRelator.xi(mx,my,q);
            s.mx = mx;
            s.my = my;
            s.q = q;
            s.phi   = 0;
            s.pNorm = 'max';
            s.print = false;
            s.iter = obj.iter;
            s.hasToCaptureImage = false;
            stressNorm = StressNormSuperEllipseComputer(s);
            obj.homog.cellVariables = stressNorm.computeCellVariables();
            obj.homog.cellVariables.mx = mx;
            obj.homog.cellVariables.my = my;
            obj.homog.cellVariables.rho = SuperEllipseParamsRelator.rho(mx,my,q);
            obj.homog.cellVariables.xi = SuperEllipseParamsRelator.xi(mx,my);
        end
        
        function obtainHomogenizerData(obj)
            obj.obtainComputedCellVariables();
            obj.obtainIntegrationVariables();
        end
        
        function obtainComputedCellVariables(obj)
            obj.obtainVolume();   
            obj.obtainCtensor();
            obj.obtainStressesAndStrain();
            obj.obtainDisplacements();
        end
        
        function obtainVolume(obj)
            v = obj.homog.cellVariables.volume;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.variables{imx,imy}.volume = v;
        end      
        
        function obtainCtensor(obj)
            Ch = obj.homog.cellVariables.Ctensor;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.variables{imx,imy}.Ctensor = Ch;
        end
        
        function obtainStressesAndStrain(obj)
            stress = obj.homog.cellVariables.tstress;            
            strain = obj.homog.cellVariables.tstrain;                        
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.variables{imx,imy}.tstress = stress;
            obj.variables{imx,imy}.tstrain = strain;            
        end
        
        function obtainDisplacements(obj)
            displ = obj.homog.cellVariables.displ;                        
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;
            obj.variables{imx,imy}.displ = displ;
        end
        
        function obtainIntegrationVariables(obj)
            intVar = obj.homog.cellVariables.integrationVar;
            imx = obj.iMxIndex;
            imy = obj.iMyIndex;            
            obj.variables{imx,imy}.integrationVar = intVar;
        end
        
    end
    
end