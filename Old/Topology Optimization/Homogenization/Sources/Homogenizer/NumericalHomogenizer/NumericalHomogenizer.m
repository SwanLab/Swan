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
        dim
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
            obj.getProblemDimensions();
            obj.createInterpolation();
            obj.createElementalDensityCreator();
            obj.obtainDensity();
            obj.createMaterialProperties()
            obj.setMaterialPropertiesInMicroProblem()
        end
        
        function buildMicroProblem(obj)
            a.fileName = obj.fileName;
            s = FemDataContainer(a);
            obj.microProblem = FEM.create(s);
        end

        function getProblemDimensions(obj)
            obj.dim = obj.microProblem.getDimensions();
        end
        
        function createInterpolation(obj)
            d = SettingsInterpolation();
            m = obj.microProblem.getMesh();
            d.interpolation = obj.interDataBase.materialInterpolation;
            d.constitutiveProperties  = obj.matDataBase.matProp;
            d.typeOfMaterial = obj.matDataBase.materialType;
            d.dim  = obj.pdim;
            d.nElem = m.nelem;
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
            d.ndim  = obj.dim.ndimf;
            d.coord = obj.microProblem.getMesh().coord;
        end
        
        function d = createFilterDataBase(obj)
            prob = obj.microProblem;
            d.shape = prob.getInterpolation().shape;
            d.conec = prob.getMesh().connec;
            d.quadr = prob.getQuadrature();
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
            obj.microProblem.solve();
            cV = obj.cellVariables;
            cV.Ch      = obj.microProblem.Chomog;
%             cV.tstress = obj.microProblem.variables.tstress;
%             cV.tstrain = obj.microProblem.variables.tstrain;
%             cV.displ   = obj.microProblem.variables.tdisp;
            obj.cellVariables = cV;
            
            
            %obj.microProblem.computeStressBasisCellProblem();
            %var = obj.microProblem.variables2printStressBasis();
        end
        
        function computeVolumeValue(obj)
            prob = obj.microProblem;
            mpMesh = prob.getMesh();
            cParams.coord  = mpMesh.coord;
            cParams.connec = mpMesh.connec;
%             mesh = Mesh_Total(cParams);
            mesh = Mesh.create(cParams);

            d = obj.volDataBase;
            s = SettingsDesignVariable();
            s.type = 'Density';
            s.mesh = mesh;%obj.microProblem.mesh;
            s.initialCase  = 'given';
            s.creatorSettings.value = obj.elemDensCr.getLevelSet();
            s.creatorSettings.ndim  = obj.dim.ndimf;
            s.creatorSettings.coord = mpMesh.coord;
            scalarPr.epsilon = 1e-3;
            scalarPr.mesh = mesh;
            s.scalarProductSettings    = scalarPr;
            d.filterParams.femSettings = d.femSettings;

            % (19/12/2023): The future idea will be to destroy
            % LevelSerCreator and use GeometricalFunction
            sLs        = s.creatorSettings;
            sLs.ndim   = s.mesh.ndim;
            sLs.coord  = s.mesh.coord;
            sLs.type   = s.initialCase;
            lsCreator  = LevelSetCreator.create(sLs);
            phi        = lsCreator.getValue();
            switch s.type
                case 'Density'
                    value = 1 - heaviside(phi);
                case 'LevelSet'
                    value = phi;
            end
            ss.fValues = value;
            ss.mesh    = s.mesh;
            ss.order   = 'P1';
            s.fun      = LagrangianFunction(ss);

            desVar = DesignVariable.create(s);
            d.filterParams.mesh = desVar.mesh;
            d.filterParams.designVarType = desVar.type;
            d.filterParams = SettingsFilter(d.filterParams);
            d.mesh = mesh;
%             d.mesh.computeMasterSlaveNodes();
            d.designVariable = desVar;

            sF            = d.filterParams.femSettings;
            sF.filterType = d.filterParams.filterType;
            sF.mesh       = d.designVariable.mesh;
            sF.test       = LagrangianFunction.create(sF.mesh, 1, 'P0');
            sF.trial      = LagrangianFunction.create(sF.mesh, 1, 'P1');
            d.femSettings.designVariableFilter = Filter.create(sF);
            d.femSettings.gradientFilter       = Filter.create(sF);

            vComputer = ShFunc_Volume(d);
            vComputer.computeFunctionFromDensity(obj.density);
            obj.cellVariables.volume = vComputer.value;
            obj.cellVariables.geometricVolume = vComputer.geometricVolume;
        end
        
        function obtainIntegrationUsedVariables(obj)
            mProb = obj.microProblem;
            intVar.geoVol = mProb.computeGeometricalVolume();
            intVar.dV     = mProb.getDvolume();
%             intVar.nstre  = obj.dim.nstre;
            intVar.ngaus  = size(intVar.dV,2);
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
                d.quad = obj.microProblem.getQuadrature();
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
            dI.mesh            = obj.microProblem.getMesh();
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

