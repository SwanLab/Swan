classdef printingVigdergauzMicroStructures < handle
    
    properties (Access = private)
        mesh
        outputName
        pdim
        postprocess
    end
    
    methods (Access = public)
        
        function init(obj)
            obj.ouputName = 'VigergauzMicroStructure';            
        end
        
        function obj = printingVigdergauzMicroStructures()
            settings = Settings('VigergauzLevelSetInput');
            
            settings.superEllipseRatio = tan(pi/6);%tan(pi/5);
            settings.volumeMicro = 0.8;
            
            
            translator = SettingsTranslator();
            translator.translate(settings);
            fileName = translator.fileName;
            settingsTopOpt = SettingsTopOptProblem(fileName);
            topOpt = TopOpt_Problem(settingsTopOpt);
            
            
            dV = topOpt.designVariable;
              
            
        end
        
    end
    
    methods (Access = private)
        
        function createPostProcess(obj)
            dB = obj.createPostProcessDataBase();
            dB.printers = obj.printers;
            postCase = 'NumericalHomogenizer';
            obj.postProcess = Postprocess(postCase,dB);
        end
        
        function dB = createPostProcessDataBase(obj)
            dI.mesh            = obj.mesh;
            dI.outName         = obj.outputName;
            dI.pdim            = obj.pdim;
            dI.ptype           = 'MICRO';
            ps = PostProcessDataBaseCreator(dI);
            dB = ps.getValue();
        end
        
    end
    
end
