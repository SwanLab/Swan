classdef NumericalHomogenizer < handle
    
    properties (Access = public)
       iter 
    end
        
    properties (SetAccess = private, GetAccess = public)
        matValues
        elemDensCr
        cellVariables
        integrationVar        
    end
    
    properties (Access = private)
        fileName
        outputName
        hasToBePrinted
        pdim
        eDensCreatType
        hasToCaptureImage = false
        
        lsDataBase
        matDataBase
        interDataBase
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
        end
        
        function compute(obj)
            obj.createMicroProblem();
            obj.computeCellVariables();
            obj.obtainIntegrationUsedVariables();
            obj.print();
            obj.captureImage();            
        end
        
        function m = getMicroProblem(obj)
            m = obj.microProblem;
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
            obj.lsDataBase     = d.levelSetDataBase;
            obj.matDataBase    = d.materialDataBase;
            obj.interDataBase  = d.materialInterpDataBase;
            obj.volDataBase    = d.volumeShFuncDataBase;
            obj.hasToCaptureImage = d.hasToCaptureImage;
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
        end        
        
        function createInterpolation(obj)
            d = SettingsInterpolation();
            d.interpolation = obj.interDataBase.materialInterpolation;
            d.constitutiveProperties  = obj.matDataBase.matProp;
            d.typeOfMaterial = obj.matDataBase.materialType;
            d.dim  = obj.pdim;
            d.nElem = obj.microProblem.mesh.nelem;
            mI  = MaterialInterpolation.create(d);
            obj.interpolation = mI;
            obj.matValues = d.constitutiveProperties;
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
            d.levelSetCreatorDataBase = dl;
            d.filterDataBase = df;
        end
        
        function d = createLevelSetCreatorDataBase(obj)
            d = obj.lsDataBase;
            d.ndim  = obj.microProblem.mesh.ndim;
            d.coord = obj.microProblem.mesh.coord;
        end
        
        function d = createFilterDataBase(obj)
            d.shape = obj.microProblem.element.interpolation_u.shape;
            d.conec = obj.microProblem.mesh.connec;
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
        
        function computeCellVariables(obj)
            obj.computeVolumeValue();
            obj.computeElasticVariables();
        end
        
        function computeElasticVariables(obj)
            obj.microProblem.computeChomog();
            cV = obj.cellVariables;
            cV.Ch      = obj.microProblem.variables.Chomog;
            cV.tstress = obj.microProblem.variables.tstress;
            cV.tstrain = obj.microProblem.variables.tstrain;
            cV.displ   = obj.microProblem.variables.tdisp; 
            obj.cellVariables = cV;
            
            
            %obj.microProblem.computeStressBasisCellProblem();
            %var = obj.microProblem.variables2printStressBasis();            
        end
        
        function computeVolumeValue(obj)
            cParams.coord  = obj.microProblem.mesh.coord;
            cParams.connec = obj.microProblem.mesh.connec;
            mesh = Mesh_Total(cParams);                        
            d = obj.volDataBase;
            s = SettingsDesignVariable();
            s.type = 'Density';            
            s.mesh = mesh;%obj.microProblem.mesh;
            s.initialCase  = 'given';
            s.creatorSettings.value = obj.elemDensCr.getLevelSet();
            s.creatorSettings.ndim  = obj.microProblem.mesh.ndim;
            s.creatorSettings.coord = obj.microProblem.mesh.coord; 
            scalarPr.epsilon = 1e-3;
            scalarPr.mesh = mesh.innerMeshOLD;
            s.scalarProductSettings    = scalarPr;
            d.filterParams.femSettings = d.femSettings;
            desVar = DesignVariable.create(s);
            d.filterParams.mesh = desVar.mesh.innerMeshOLD;
            d.filterParams.designVarType = desVar.type;
            d.filterParams = SettingsFilter(d.filterParams);
            d.mesh = mesh.innerMeshOLD;
            d.mesh.computeMasterSlaveNodes();
            d.designVariable = desVar;
            vComputer = ShFunc_Volume(d);
            vComputer.computeFunctionFromDensity(obj.density);
            obj.cellVariables.volume = vComputer.value;
            obj.cellVariables.geometricVolume = vComputer.geometricVolume;
        end
        
        function mesh = setMasterSlaveNodes(obj,mesh)
                       
        end
               
        function obtainIntegrationUsedVariables(obj)        
           intVar.nstre  = obj.microProblem.element.getNstre();
           intVar.geoVol = obj.microProblem.computeGeometricalVolume();
           intVar.ngaus  = obj.microProblem.element.quadrature.ngaus;
           intVar.dV     = obj.microProblem.geometry.dvolu;
           obj.integrationVar = intVar;
        end
        
        function print(obj)
            if obj.hasToBePrinted
                obj.createPrintersNames();
                obj.createPostProcess();
                d.var2print = obj.elemDensCr.getFieldsToPrint;
                d.var2print{end+1} = obj.microProblem;
               % obj.microProblem.variables.var2print = obj.microProblem.variablesStressBasis.var2print;
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
            obj.printers{end+1} = 'HomogenizedTensorStressBasis';
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
            dI.pdim            = obj.pdim;
            dI.ptype           = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
        end
        
        function captureImage(obj)
            if obj.hasToCaptureImage
                i = obj.iter;
                f = obj.resFile;
                outPutNameWithIter = [obj.outputName,num2str(i)];
                inputFileName = fullfile('Output',f,[f,num2str(i),'.flavia.res']);
                s.fileName = f;
                s.outPutImageName = outPutNameWithIter;
                s.inputFileName = inputFileName;
                imageCapturer = GiDImageCapturer(s);
                imageCapturer.capture();
            end
        end
        
    end
    
end

