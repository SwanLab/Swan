classdef Postprocess < handle
       
    properties (Access = protected)                  
        outFileName
        
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
               
        function  print(obj,iter,results)
            obj.printMeshFile(iter);            
            obj.printResFile(iter,results)
        end
        
        function r = getResFile(obj)
            r = obj.resultPrinter.getFieldName();
        end
        
        function setIter(obj,i)
            obj.resultPrinter.setIter(i);
        end
        
    end
       
    methods (Access = private)
        
        function init(obj,d)
            obj.outFileName = d.outFileName;            
            obj.mshDataBase = obj.computeDataBaseForMeshFile(d);                       
            obj.resDataBase = obj.computeDataBaseForResFile(d); 
        end
        
        function createMeshPrinter(obj)
            obj.meshPrinter = MeshPrinter();            
        end
        
        function createResultPrinter(obj,postCase)
            factory  = ResultsPrinterFactory();
            obj.resultPrinter = factory.create(postCase,obj.resDataBase);
        end
        
        function printResFile(obj,iter,results)           
            obj.resultPrinter.print(iter,results);            
        end
        
       	function printMeshFile(obj,iter)
            d = obj.mshDataBase;
            d.iter = iter;
            obj.meshPrinter.print(d);
        end
        
        function createOutputDirectory(obj)
            path = pwd;
            dir = fullfile(path,'Output',obj.outFileName);
            if ~exist(dir,'dir')
                mkdir(dir)
            end            
        end
    end
    
     methods (Access = private, Static)
        
        function d = computeDataBaseForResFile(dI)
            d.testName = dI.outFileName;
            d.etype = dI.etype;
            d.ptype = dI.ptype;
            if isfield(dI,'ngaus')
                d.ngaus = dI.ngaus;
            else
                d.ngaus = [];
            end
            d.ndim = dI.ndim;
            if isfield(dI,'posgp')
                d.posgp = dI.posgp;
            else
                d.posgp = [];
            end
            
            if isfield(dI,'ShapeNames')
                d.ShapeNames = dI.ShapeNames;
            end
            
            if isfield(dI,'optimizer')
                d.optimizer = dI.optimizer;
            end
            
            if isfield(dI,'printMode')
                d.printMode = dI.printMode;                
            end
            
            if isfield(dI,'hasGaussData')
                d.hasGaussData = dI.hasGaussData;                
            end   
            d.gaussDescriptor = 'Guass up?';
        end
        
        function d = computeDataBaseForMeshFile(dI)
            d.coordinates = dI.coordinates;
            d.connectivities = dI.connectivities;
            d.testName = dI.outFileName;
            d.npnod = dI.npnod;
            d.pdim = dI.pdim;
            d.nnode = dI.nnode;
            d.nelem = dI.nelem;
            d.ndim = dI.ndim;
            d.etype = dI.etype;
        end
 
    end

end
