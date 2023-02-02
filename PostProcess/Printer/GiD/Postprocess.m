classdef Postprocess < handle
       
    properties (Access = protected)
        outFileName
        resultsDir
        
        resultPrinter
        meshPrinter
        
        mshDataBase
        resDataBase
    end
    
    
     methods (Access = public)
 
        function obj = Postprocess(postCase,d)
            obj.init(d);
            obj.createOutputDirectory();
            obj.createMeshPrinter();
            obj.createResultPrinter(postCase);
        end
        
        function print(obj,iter,d)
            obj.printMeshFile(iter);
            obj.printResFile(iter,d)
        end
        
        function r = getResFile(obj)
            r = obj.outFileName;
        end
        
    end
       
    methods (Access = private)
        
        function init(obj,d)
            obj.outFileName = d.outFileName;
            obj.createResultsDirName();
            obj.computeDataBaseForMeshFile(d);
            obj.computeDataBaseForResFile(d);
        end
        
        function createMeshPrinter(obj)
            obj.meshPrinter = MeshPrinter();
        end
        
        function createResultPrinter(obj,postCase)
            obj.resultPrinter  = ResultsPrinter.create(postCase,obj.resDataBase);
        end
        
        function printResFile(obj,iter,d)
            obj.resultPrinter.print(iter,d);
        end
        
       	function printMeshFile(obj,iter)
            d = obj.mshDataBase;
            d.iter = iter;
            obj.meshPrinter.print(d);
        end
        
        function createOutputDirectory(obj)
            dir = obj.resultsDir;
            if ~exist(dir,'dir')
                mkdir(dir)
            end
        end
        
        function computeDataBaseForResFile(obj,dI)
            obj.resDataBase = dI;
            obj.resDataBase.resultsDir = obj.resultsDir;
        end
        
        function computeDataBaseForMeshFile(obj,dI)
            d.coordinates    = dI.coordinates;
            d.connectivities = dI.connectivities;
            d.outFileName    = dI.outFileName;
            d.npnod          = dI.npnod;
            d.pdim           = dI.pdim;
            d.nnode          = dI.nnode;
            d.nelem          = dI.nelem;
            d.ndim           = dI.ndim;
            d.etype          = dI.etype;
            d.resultsDir     = obj.resultsDir;
            obj.mshDataBase = d;
        end
        
        function createResultsDirName(obj)
           path = pwd;
           obj.resultsDir = fullfile(path,'Output',obj.outFileName);
        end
         
     end

end
