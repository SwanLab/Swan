classdef NumericalHomogenizer < handle 
    

    properties (Access = private)
        fileName
        fileNameWithPath
        hasToBePrinted
        outputName  
        matProp
        matValues
        Ptensor
        volume
        
        densityPrinter
        
        interpolation
        densityPostProcess
        
        iter
    end
    
    properties (Access = protected)
        microProblem
        density
        Ch
        setting
    end
    
    methods (Access = public)
        
        function p = getAmplificatorTensor(obj)
            p = obj.Ptensor;
        end
        
        function c = getCh(obj)
            c = obj.Ch;        
        end
        
        function m = getMaterialValues(obj)
            m = obj.matValues;
        end
        
        function v = getVolume(obj)
            v = obj.volume;
        end
        
        
    end
    
    methods (Access = protected)

        function init(obj,outputName,print,iter)
            obj.outputName = outputName;
            obj.hasToBePrinted = print;
            obj.iter = iter;
            obj.loadFileName();
            obj.loadTestData();
            obj.loadDimension();
        end
        
        function generateMicroProblem(obj)
            obj.createInterpolation();
            obj.createMicroProblem();
            obj.createDensity()
            obj.createMaterialProperties()
            obj.setMaterialPropertiesInMicroProblem()
        end             
        
        function createDensityPrinter(obj)
            quad = obj.microProblem.element.quadrature;
            mesh = obj.microProblem.mesh;
            obj.densityPrinter = DensityPrinter(quad,mesh);
        end        
        
        function print(obj)
            if obj.hasToBePrinted
                d = obj.density;            
                outn = obj.outputName;
                it   = obj.iter;                
                obj.densityPrinter.print(d,outn,it)
            end
        end         
        
        function computeHomogenizedVariables(obj)
            obj.computeVolumeValue()
            obj.computeHomogenizedTensor()
            obj.computeAmplificator()
        end        
        
    end
    
    methods (Access = private)
        
        
        function loadFileName(obj)
            obj.fileName = 'test_microFineFine';
        end
        
        function loadTestData(obj)
            obj.fileNameWithPath = strcat('./Input/',obj.fileName);
            obj.setting = Settings(obj.fileNameWithPath);
        end
        
        function loadDimension(obj)
            obj.setting.pdim = '2D';
        end
                
        function createInterpolation(obj)
            matVal          = obj.setting.TOL;
            material        = obj.setting.material;
            InterpFunction  = obj.setting.method;
            dim             = obj.setting.pdim;
            interp          = Material_Interpolation.create(matVal,material,InterpFunction,dim);
            obj.interpolation = interp;
            obj.matValues = matVal;
        end
        
        function createMicroProblem(obj)
            obj.microProblem = Elastic_Problem_Micro(obj.setting.filename);
            obj.microProblem.preProcess();
        end
        
        function createMaterialProperties(obj)
            d = obj.density;
            obj.matProp= obj.interpolation.computeMatProp(d);
        end
        
        function setMaterialPropertiesInMicroProblem(obj)
            obj.microProblem.setMatProps(obj.matProp);
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
        
        function computeVolumeValue(obj)
            vComputer = ShFunc_Volume(obj.setting);
            dens = obj.density;
            vol = vComputer.computeCost(dens);
            obj.volume = vol;
        end        
        
    end
    
    methods (Access = protected, Abstract)
        createDensity(obj)
    end
    
end

