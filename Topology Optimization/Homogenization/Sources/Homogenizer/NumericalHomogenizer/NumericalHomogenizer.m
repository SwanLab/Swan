classdef NumericalHomogenizer < handle
    
    
    properties (SetAccess = private, GetAccess = public)
        Ptensor
        Ch
        volume
        matValues
        elemDensCr
    end
    
    properties (Access = private)
        fileName
        outputName
        hasToBePrinted
        iter
        pdim
        eDensCreatType
        hasToCaptureImage = false
        
        lsDataBase
        interpDataBase
        volDataBase
        
        matProp
        
        postProcess
        
        microProblem
        density
        levelSet
        resFile
        
        printers
        interpolation
    end
    
    methods (Access = public)
        
        function obj = NumericalHomogenizer(d)
            obj.init(d);
            obj.createMicroProblem();
            obj.computeHomogenizedVariables();
            obj.print();
            obj.captureImage();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,d)
            obj.fileName       = d.testName;
            obj.outputName     = d.outFileName;
            obj.hasToBePrinted = d.print;
            obj.iter           = d.iter;
            obj.pdim           = d.pdim;
            obj.eDensCreatType = d.elementDensityCreatorType;
            obj.lsDataBase     = d.levelSetCreatorParams;
            obj.interpDataBase = d.interpParams;
            obj.volDataBase    = d.volumeShFuncParams;
        end
        
        function createMicroProblem(obj)
            obj.buildMicroProblem();            
            obj.createInterpolation();
            obj.createElementalDensityCreator();
            obj.obtainDensity();
            obj.createMaterialProperties()
            obj.setMaterialPropertiesInMicroProblem()
        end
        
        function buildMicroProblem(obj)
            obj.microProblem = Elastic_Problem_Micro(obj.fileName);
            obj.microProblem.preProcess();
        end        
        
        function createInterpolation(obj)
            mI  = Material_Interpolation.create(obj.interpDataBase);
            obj.interpolation = mI;
        end
        
        function createElementalDensityCreator(obj)
            type = obj.eDensCreatType;
            de   = obj.createElementalDensityCreatorDataBase();
            edc  = ElementalDensityCreator.create(type,de);
            obj.elemDensCr = edc;
        end
        
        function d = createElementalDensityCreatorDataBase(obj)
            dl = obj.createLevelSetCreatorDataBase();
            df = obj.createFilterDataBase();
            d = SettingsElementalDensity();
            d.levelSetCreatorParams = dl;
            d.filterParams = df;
        end
        
        function d = createLevelSetCreatorDataBase(obj)
            d = obj.lsDataBase;
            d.ndim  = obj.microProblem.mesh.ndim;
            d.coord = obj.microProblem.mesh.coord;
        end
        
        function d = createFilterDataBase(obj)
            d = SettingsFilterP0();
            d.shape = obj.microProblem.element.interpolation_u.shape;
            d.conec = obj.microProblem.geometry.interpolation.T;
            d.quadr = obj.microProblem.element.quadrature;
        end        
        
        function obtainDensity(obj)
            obj.density = obj.elemDensCr.getDensity();
        end
        
        function createMaterialProperties(obj)
            d = obj.density;
            obj.matProp = obj.interpolation.computeMatProp(d);
        end
        
        function setMaterialPropertiesInMicroProblem(obj)
            obj.microProblem.setMatProps(obj.matProp);
        end
        
        function computeHomogenizedVariables(obj)
            obj.computeVolumeValue();
            obj.computeHomogenizedTensor();
            obj.computeAmplificator();
            obj.computeGeneralizedAmplificatorTensor();
        end
        
        function computeHomogenizedTensor(obj)
            obj.microProblem.computeChomog();
            obj.Ch = obj.microProblem.variables.Chomog;
        end
        
        function computeAmplificator(obj)
            Pv = obj.microProblem.computeAmplificator();
            P = SymmetricFourthOrderPlaneStressVoigtTensor();
            P.setValue(Pv);
            obj.Ptensor = P;
        end
        
        function computeGeneralizedAmplificatorTensor(obj)
            Pg = obj.microProblem.computeGeneralizedAmplificator();
        end
        
        function computeVolumeValue(obj)
            d = obj.volDataBase;
            vComputer = ShFunc_Volume(d);
            dens = obj.density;
            vol = vComputer.computeCost(dens);
            obj.volume = vol;
        end
        
        function print(obj)
            if obj.hasToBePrinted
                obj.createPrintersNames();
                obj.createPostProcess();
                d.var2print = obj.elemDensCr.getFieldsToPrint;
                d.var2print{end+1} = obj.microProblem;
                d.quad = obj.microProblem.element.quadrature;
                obj.postProcess.print(obj.iter,d);
                obj.resFile = obj.postProcess.getResFile();
            end
        end
        
        function createPrintersNames(obj)
            type = obj.eDensCreatType;
            f = ElementalDensityCreatorFactory();
            obj.printers = f.createPrinters(type);
            obj.printers{end+1} = 'HomogenizedTensor';
        end
        
        function createPostProcess(obj)
            dB = obj.createPostProcessDataBase();
            dB.printers = obj.printers;
            postCase = 'NumericalHomogenizer';
            obj.postProcess = Postprocess(postCase,dB);
        end
        
        function dB = createPostProcessDataBase(obj)
            dI.mesh            = obj.microProblem.mesh;
            dI.outName         = obj.outputName;
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
        end
        
        function captureImage(obj)
            if obj.hasToCaptureImage
                i = obj.iter;
                f = obj.resFile;
                outPutNameWithIter = [obj.outputName,num2str(i)];
                inputFileName = fullfile('Output',f,[f,num2str(i),'.flavia.res']);
                GiDImageCapturer(f,outPutNameWithIter,inputFileName);
            end
        end
        
    end
    
end

