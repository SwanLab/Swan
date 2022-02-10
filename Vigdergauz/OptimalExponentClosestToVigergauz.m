classdef OptimalExponentClosestToVigergauz < handle
    
    properties (Access = public)
        qOpt
        error
        mx
        my
        frames
        qIter
    end
    
    properties (Access = private)
        txi
        rho
        comparator
        problem
        savingFrames
    end
    
    methods (Access = public)
        
        function obj = OptimalExponentClosestToVigergauz(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            p.objective = @(q) obj.vigdergauzSuperEllipseDistance(q);
            p.x1 = 2;
            p.x2 = 32;
            p.solver = 'fminbnd';
            p.options = optimset('Display','iter','TolX',1e-8,'MaxIter',1000);
            obj.problem = p;
            [x,fsol] = fminbnd(obj.problem);
            obj.qOpt = x;
            obj.error = fsol;
            obj.mx = obj.comparator.mx;
            obj.my = obj.comparator.my;
        end
        
    end
    
   methods (Access = private)
       
       function init(obj,cParams)
          obj.rho = cParams.rho;
          obj.txi = cParams.txi;
          obj.savingFrames = cParams.savingFrames;
          obj.comparator = VigdergauzSuperEllipseComparator;
       end
       
        function d = vigdergauzSuperEllipseDistance(obj,q)
            obj.comparator.computeComparison(obj.txi,obj.rho,q);
            if obj.savingFrames
                obj.comparator.plotDifference();
                obj.frames{end+1} = obj.comparator.frame;
                obj.qIter{end+1} = q;
            end
            d = obj.comparator.rhoDifferenceNorm;
        end
       
   end
    
end