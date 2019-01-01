classdef testFEMPrinting < ...
            testNotShowingError & testFemComputation
    
    properties (Access = protected)
        testName = 'test2d_quad';  
        postProcessor = 'Elasticity';
        filesHaveChanged
        fileOutputName
    end
    
    properties (Access = private)
       dataBase
       iter = 0;
    end
    
    methods (Access = public)
        
        function obj = testFEMPrinting()
            obj.init()
            obj.print()
            obj.compareFiles()
        end
        
    end
    
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            hasPassed = ~obj.filesHaveChanged;
        end
        
        function selectComputedVar(obj)
        end
        
    end

    methods (Access = private)
        
        function init(obj)
            obj.fileOutputName = ['testFemPrinting'];
        end
        
        function compareFiles(obj)
            hasMshChanged = obj.compareFile('.msh');
            hasResChanged = obj.compareFile('.res');
            obj.filesHaveChanged = hasMshChanged || hasResChanged;
        end
        
        function hasChanged = compareFile(obj,extension)
            out   = obj.fileOutputName;
            fullOutName = [out,num2str(obj.iter),'.flavia',extension];
            savedPrintedFile = fullfile('tests','PrintingTests','PrintedFiles',fullOutName);
            outputFile = fullfile('Output',out,fullOutName);
            command = ['diff ', savedPrintedFile, ' ', outputFile];
            [hasChanged,~] = system(command);
        end       
                
        function print(obj)            
            obj.createPostProcessDataBase();
            postprocess = Postprocess(obj.postProcessor);                                            
            postprocess.print(obj.dataBase);
        end
        
        function createPostProcessDataBase(obj)
            d.mesh    = obj.fem.mesh;
            d.fields  = obj.fem.variables;
            d.outName = obj.fileOutputName;
            d.quad    = obj.fem.element.quadrature;
            d.iter    = 0;
            hasGaussInfo = true;
            ps = PostProcessDataBaseCreator.create(hasGaussInfo,d);
            obj.dataBase = ps.getValue();   
        end

    end
    
end

