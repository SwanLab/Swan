classdef RegularizedPerimeterComputer < handle
    
    properties (Access = public)
        perimeters
        perimetersGradient
        regularizedDensity
        perimeterShapeFunctions
        epsilons
    end
    
    properties (Access = private)
        inputFile
        backgroundMesh
        designVariable
        scale
        nEpsilon
        epsilon
        printing
        plotting
        capturingImage
        perimeterShapeFunction
        outputFigureName
        perimeterType
        isRobinTermAdded
    end
    
    methods (Access = public)
        
        function obj = RegularizedPerimeterComputer(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.createEpsilonValues();
            obj.computePerimeters();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.inputFile        = cParams.inputFile;
            obj.backgroundMesh   = cParams.backgroundMesh;
            obj.scale            = cParams.scale;
            obj.designVariable   = cParams.designVariable;
            obj.outputFigureName = cParams.outputFigureName;
            obj.plotting         = cParams.plotting;
            obj.printing         = cParams.printing;
            obj.capturingImage   = cParams.capturingImage;
            obj.perimeterType    = cParams.perimeterType;
            obj.isRobinTermAdded = cParams.isRobinTermAdded;
        end
        
        function createEpsilonValues(obj)
            epsmin = obj.backgroundMesh.computeMeanCellSize;
            epsmax = obj.backgroundMesh.computeCharacteristicLength();
            obj.nEpsilon = min(6,ceil(log2(epsmax/epsmin)));
            obj.epsilons = epsmin*(2.^((1:obj.nEpsilon) - 1));
        end
        
        function computePerimeters(obj)
            obj.perimeters = zeros(obj.nEpsilon,1);
            for iepsilon = 1:obj.nEpsilon
                obj.epsilon = obj.epsilons(iepsilon);
                obj.computeRegularizedPerimeter(iepsilon);
                obj.plotDensity();
                obj.printDensity(iepsilon);
                obj.captureImage(iepsilon);
            end
        end
        
        function computeRegularizedPerimeter(obj,iepsilon)
            obj.createPerimeterShapeFunction();
            obj.perimeterShapeFunction.computeFunctionAndGradient();
            per    = obj.perimeterShapeFunction.value;
            dPer   = obj.perimeterShapeFunction.gradient;
            rhoReg = obj.perimeterShapeFunction.regularizedDensity;
            obj.perimeters(iepsilon) = per*obj.perimeterShapeFunction.value0;
            obj.perimetersGradient(:,iepsilon) = dPer;
            obj.regularizedDensity(:,iepsilon) = rhoReg;
            obj.perimeterShapeFunctions{iepsilon} = obj.perimeterShapeFunction;
        end
        
        function createPerimeterShapeFunction(obj)
            s = obj.createPerimeterParams();
            shFunc = ShFunc_Perimeter(s);
            obj.perimeterShapeFunction = shFunc;
        end
        
        function s = createPerimeterParams(obj)
            sC.inputFile        = obj.inputFile;
            sC.mesh             = obj.backgroundMesh;
          %  sC.designVariable.value      = obj.designVariable;
          %  sC.designVariable.nVariables = 1;
            sC.designVariable   = obj.createDesignVariable();
            sC.epsilon          = obj.epsilon;
            sC.scale            = obj.scale;
            sC.type             = obj.perimeterType;
            sC.isRobinTermAdded = obj.isRobinTermAdded;
            fCreator = PerimeterParamsCreator(sC);
            s = fCreator.perimeterParams;
        end
        
        function dV = createDesignVariable(obj)
          %  s.scalarProductSettings.femSettings = [];
          %  s.epsilon = [];  
          %  ss.type = 'full';
          %  sLs = SettingsLevelSetCreator;
          %  s = sLs.create(ss); 
            s.mesh      = obj.backgroundMesh;
            s.inputFile = obj.inputFile;
            s.scale     = obj.scale;
            s.type      = 'LevelSet';
            d = DesignVariableCreatorSettings(s);
            s = d.create();
            dV = DesignVariable.create(s);                                    
            dV.update(obj.designVariable);
        end
        
        function plotDensity(obj)
            s.mesh      = obj.backgroundMesh;
            s.inputFile = obj.inputFile;
            s.scale     = obj.scale;
            s.density   = obj.perimeterShapeFunction.regularizedDensity;
            plotter = DensityPlotterForPerimeter(s);
            if obj.plotting
                plotter.plot();
            end
        end
        
        function printDensity(obj,iepsilon)
            s.inputFile = obj.inputFile;
            s.iter      = iepsilon;
            s.mesh      = obj.backgroundMesh;
            s.perimeter   = obj.perimeterShapeFunction;
            printer = DensityPrinterForPerimeter(s);
            printer.print();
        end
        
        function captureImage(obj,iepsilon)
            if obj.capturingImage
                i = iepsilon;
                f = obj.inputFile;
                outPutNameWithIter = [obj.outputFigureName,'Epsilon',num2str(i)];
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