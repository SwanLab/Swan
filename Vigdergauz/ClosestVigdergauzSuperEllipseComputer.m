classdef ClosestVigdergauzSuperEllipseComputer < handle
    
    properties (Access = public)
        xopt
        error
    end
    
    properties (Access = private)
        txi
        rho
        q
        comparator
        problem
        frames
    end
    
    methods (Access = public)
        
        function  obj = ClosestVigdergauzSuperEllipseComputer()
            obj.init();
            obj.solveProblem();
            obj.writeOptimizationVideo();
            obj.printVigdergauzSuperEllipse()
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.comparator = VigdergauzSuperEllipseComparator;
            obj.txi = pi/8;
            obj.rho = 0.8;
        end
        
        function solveProblem(obj)
            p.objective = @(q) obj.vigdergauzSuperEllipseDistance(q);
            p.x1 = 2;
            p.x2 = 32;
            p.solver = 'fminbnd';
            p.options = optimset('Display','iter','TolX',1e-8,'MaxIter',1000);
            obj.problem = p;
            [x,fsol] = fminbnd(obj.problem);
            obj.xopt = x;
            obj.q = x;
            obj.error = fsol;
        end
        
        function d = vigdergauzSuperEllipseDistance(obj,q)
            obj.comparator.computeComparison(obj.txi,obj.rho,q);
            obj.frames{end+1} = obj.comparator.frame;
            d = obj.comparator.rhoDifferenceNorm;
        end
        
        function writeOptimizationVideo(obj)
            v = VideoWriter('SuperEllipseToVigdergauz.avi');
            open(v);
            nf = 6;
            for iframe = 1:nf*(numel(obj.frames)-2)
                i = floor(iframe/nf)+1;
                frame = obj.frames{i};
                writeVideo(v,frame);
            end            
            close(v);
        end
        
        function printVigdergauzSuperEllipse(obj)
            s.mx = obj.comparator.mx;
            s.my = obj.comparator.my;
            s.txi = obj.txi;
            s.rho = obj.rho;
            s.q   = obj.q;
            stressPrinter = VigdergauzSuperEllipseStressPrinter(s);
            stressPrinter.print();            
        end
        
    end
    
end

