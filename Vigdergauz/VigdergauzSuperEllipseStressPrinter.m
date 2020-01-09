classdef VigdergauzSuperEllipseStressPrinter < handle
       
   properties (Access = private)
       mx
       my
       txi
       rho
       q
       homogDataBase
       homog
       fileCase
   end
    
   methods (Access = public)
        
       function obj = VigdergauzSuperEllipseStressPrinter(cParams)
           obj.init(cParams)
       end       
       
       function print(obj)
           obj.printStressSuperEllipse();
           obj.printStressVigdergauz();           
       end
       
   end
   
   methods (Access = private)
       
       function init(obj,cParams)
           obj.txi = cParams.txi;
           obj.rho = cParams.rho;
           obj.q   = cParams.q;
           obj.mx = cParams.mx;
           obj.my = cParams.my;
       end
       
       function printStressSuperEllipse(obj)
           obj.createNumericalHomogenizerDataBaseForSuperEllipse();
           obj.createNumericalHomogenizer();
       end
       
       function printStressVigdergauz(obj)
           obj.createNumericalHomogenizerDataBaseForVigdergauz();
           obj.createNumericalHomogenizer();
       end       
       
        function createNumericalHomogenizerDataBaseForSuperEllipse(obj)
            obj.fileCase = 'RVE_Square_Triangle_FineFine';
            defaultDB = NumericalHomogenizerDataBase([obj.fileCase,'.m']);
            dB = defaultDB.dataBase;
            dB.print                         = true;
            dB.hasToCaptureImage = false;            
            dB.outFileName                   = 'SuperEllipseLevelSetInput';            
            dB.levelSetDataBase.vigdergauzDataBase.volumeMicro = obj.rho;                                    
            dB.levelSetDataBase.type = 'smoothRectangle';
            dB.levelSetDataBase.widthH = obj.mx;
            dB.levelSetDataBase.widthV = obj.my;
            dB.levelSetDataBase.pnorm = obj.q;
            obj.homogDataBase = dB;
        end  
       
        function createNumericalHomogenizerDataBaseForVigdergauz(obj)
            obj.fileCase = 'RVE_Square_Triangle_FineFine';
            defaultDB = NumericalHomogenizerDataBase([obj.fileCase,'.m']);
            dB = defaultDB.dataBase;
            dB.print                         = true;
            dB.hasToCaptureImage = false;            
            dB.outFileName                   = 'VigergauzLevelSetInput';            
            dB.levelSetDataBase.vigdergauzDataBase.volumeMicro = obj.rho;                                    
            dB.levelSetDataBase.type = 'Vigdergauz';
            dB.levelSetDataBase.vigdergauzDataBase.superEllipseRatio = tan(obj.txi);
            dB.levelSetDataBase.vigdergauzDataBase.type = 'VolumeAndRatio';   
            obj.homogDataBase = dB;
        end  
        
        function createNumericalHomogenizer(obj)
            d = obj.homogDataBase;
            obj.homog = NumericalHomogenizer(d);
            obj.homog.iter = 0;
            obj.homog.compute();                                
        end
       
       
   end
       
end