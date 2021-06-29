classdef  SmoothRectangleVademecumComputerForGivenVolume < ...
        VademecumComputerForGivenVolume
    
    properties (Access = protected)
       qValue 
    end
    
    methods (Access = public)
        
        function obj = SmoothRectangleVademecumComputerForGivenVolume(d)
            obj.init(d);
            obj.fileName = 'SmoothRectangle';
            obj.qValue = obj.qNorm;
        end
        
    end
    
    methods (Access = protected)        
        
        function findInclusionLengthForCertainVolume(obj)
            obj.print = false;
            x0 = [obj.computeInclusionLengthForRectangle,0.99];
            problem.objective = @(x) obj.fzero(x);
            problem.x0 = x0;
            problem.solver = 'fzero';
            problem.options = optimset('Display','iter','TolX',1e-5);
            [xroot,~] = fzero(problem);
            obj.my = xroot;
        end      
        
        function f = fzero(obj,x)
            obj.computeCellVariables(x);
            vol = obj.obtainVolume();
            f = (vol - obj.volume);
        end           
        
        function vol = obtainVolume(obj)
            fileName = [obj.prefixName,obj.fileName];
            matFile = [fileName,'.mat'];
            file2load = fullfile('Output',fileName,matFile);
            a  = load(file2load);
            vad = a.d;   
            vol = vad.variables{1,1}.volume();            
        end
        
    end
    
    
    
end