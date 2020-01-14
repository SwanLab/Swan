classdef RegularizedPerimeterComputer < handle
    
    properties (Access = public)
        perimeters
        epsilons
    end
    
    properties (Access = private)
        inputFile
        mesh
        designVariable
        scale
        nEpsilon
        epsilon
        printing
        plotting
        capturingImage
        perimeterShapeFunction
        outputFigureName
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
            obj.mesh             = cParams.mesh;
            obj.scale            = cParams.scale;
            obj.designVariable   = cParams.designVariable;
            obj.outputFigureName = cParams.outputFigureName;
            obj.plotting         = cParams.plotting;
            obj.printing         = cParams.printing;
            obj.capturingImage   = cParams.capturingImage;
        end
        
        function createEpsilonValues(obj)
            m = obj.designVariable.mesh;
            epsmin = m.computeMeanCellSize;
            epsmax = m.innerMeshOLD.computeCharacteristicLength();
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
            obj.perimeterShapeFunction.computeCostAndGradient();
            per = obj.perimeterShapeFunction.value;
            obj.perimeters(iepsilon) = per;
        end
        
        function createPerimeterShapeFunction(obj)
            s = obj.createPerimeterParams();
            shFunc = ShFunc_Perimeter(s);
            obj.perimeterShapeFunction = shFunc;
        end
        
        function s = createPerimeterParams(obj)
            sC.inputFile      = obj.inputFile;
            sC.mesh           = obj.mesh;
            sC.designVariable = obj.designVariable;
            sC.epsilon        = obj.epsilon;
            sC.scale          = obj.scale;
            fCreator = PerimeterParamsCreator(sC);
            s = fCreator.perimeterParams;
        end
        
        function plotDensity(obj)
            s.mesh      = obj.mesh;
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
            s.mesh      = obj.designVariable.mesh;
            s.density   = obj.perimeterShapeFunction.regularizedDensity;
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