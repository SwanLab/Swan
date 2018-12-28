classdef Printer < handle
    
    properties (Access = private)
       resultsFileName
    end
    
    
    properties (Access = protected)
       postProcess
    end
    
    methods (Access = public)
        
        function print(obj,variable,outname,iter,quad,mesh)
            obj.postProcess = Postprocess(postCase);
            dI.x = variable;
            dI.fileOutputName = outname;
            dI.iter = iter;
            dI.mesh = mesh;
            dI.quad = quad;
            d = obj.createPostProcessDataBaseStructre(dI);
            obj.postProcess.print(d)
            obj.resultsFileName = obj.PostProcess.getResFile();
        end
        
        function r = getResFile(obj)
            r = obj.resultsFileName;
        end
        
    end
            
    
    methods (Access = private)
        
        function d = createPostProcessDataBaseStructre(obj,dI)
       
        end
        
    end
    
    methods (Abstract,Access =  protected, Static)
     createPostProcess(obj)        
    end
    
end
