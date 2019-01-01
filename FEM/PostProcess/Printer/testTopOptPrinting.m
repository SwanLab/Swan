classdef testTopOptPrinting < testNotShowingError ...
                                  & testTopOptComputation ...
                                  & testPrintingDescriptor
    
    properties (Access = protected)
      filesHaveChanged
      iter
      fields
    end
    
    properties (Access = protected, Abstract)
       postProcessor 
    end
    
    methods (Access = protected)
        
        function obj = testTopOptPrinting()
            obj.init();
            obj.print();
            obj.compareFiles();
        end
        
        function init(obj)
            obj.iter = 0;
        end
        
        function print(obj)
            obj.computeFields()
            postprocess = Postprocess(obj.postProcessor);
            d = obj.createPostProcessDataBaseStructre();                     
            postprocess.print(d);
        end
        
        function d = createPostProcessDataBaseStructre(obj)
            dI.mesh    = obj.topOpt.mesh;
            dI.fields  = obj.fields;
            dI.outName = obj.fileOutputName;
            dI.iter    = obj.iter;
            hasGaussData = false;
            ps = PostProcessDataBaseCreator.create(hasGaussData,dI);
            d = ps.getValue();              
        end
        
    end
    
    methods (Access = protected, Abstract)
       computeFields(obj) 
    end
    
end

