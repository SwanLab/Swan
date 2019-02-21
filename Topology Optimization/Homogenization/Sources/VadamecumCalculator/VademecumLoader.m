classdef VademecumLoader < handle
    
    properties (SetAccess = private, GetAccess = public)
        dataBase        
    end
    
    properties (Access = private)
        index
        path
        iIndex        
        jIndex   
        tensor        
        microName
        tensorCase  
        tensName        
    end
    
    methods (Access = public)
        
        function obj = VademecumLoader(d)
            obj.init(d)
            obj.loadMxMy();
            obj.loadHomogenizedTensor()
            obj.loadAmplificatorTensor()
        end
        
    end    
    
    methods (Access = private)
        
        function init(obj,d)
            obj.index = d.index;
            obj.path = d.path;
            obj.microName = d.microName;
        end
        
        function loadMxMy(obj)
            m1Path = fullfile(obj.path,'m1m2.txt');
            obj.dataBase.mx = load(m1Path);
            obj.dataBase.my = obj.dataBase.mx;
            obj.dataBase.nmx = length(obj.dataBase.mx);
            obj.dataBase.nmy = length(obj.dataBase.my);
        end
        
        function loadHomogenizedTensor(obj)
            obj.tensorCase = 'C_';
            obj.tensName = 'ConstitutiveTensor';
            obj.loadTensors();
            obj.dataBase.C = obj.tensor;
        end
        
        function loadAmplificatorTensor(obj)
            obj.tensorCase = 'P_';
            obj.tensName = 'AmplificatorTensor';
            obj.loadTensors();
            obj.dataBase.P = obj.tensor;
        end
        
        function loadTensors(obj)
            nPlot = length(obj.index);
            for iplot = 1:nPlot
                obj.obtainIJindex(iplot);
                obj.loadTensor();
            end
        end
        
       function obtainIJindex(obj,iplot)
            obj.iIndex = obj.index(iplot,1);
            obj.jIndex = obj.index(iplot,2);
       end              
        
        function loadTensor(obj)
            mN    = obj.microName;
            tName = obj.tensName;
            i = obj.iIndex;
            j = obj.jIndex;
            name = [obj.tensorCase,num2str(i),num2str(j)];
            fullPath = fullfile(obj.path,mN,tName,[mN,name,'.txt']);
            obj.tensor(i,j,:,:) = load(fullPath);
        end
                
    end
    
end