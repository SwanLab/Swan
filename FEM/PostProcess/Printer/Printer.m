classdef Printer < handle
    
    properties (Access = private)
       mesh 
       resultsFileName
    end
    
    
    properties (Access = protected)
       quadrature        
       PostProcess
       iter
    end
    
    methods (Access = public)
        
        function print(obj,variable,outname,iter)
            m = obj.mesh;
            obj.PostProcess.print(m,variable,iter,outname)
            obj.resultsFileName = obj.PostProcess.getResFile();
        end
        
        function r = getResFile(obj)
            r = obj.resultsFileName;
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,quad,mesh)
            obj.quadrature = quad;
            obj.mesh = mesh;
            obj.createPostProcess()
        end
        
    end
    
    methods (Abstract,Access =  protected, Static)
     createPostProcess(obj)        
    end
    
end
