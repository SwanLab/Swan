classdef ClosestVigdergauzSuperEllipseComputer < handle
    
    properties (Access = public)
        xopt
        error
    end
    
    properties (Access = private)
        txi
        rho
        q
        problem
        frames
        optimalExponent
    end
    
    methods (Access = public)
        
        function  obj = ClosestVigdergauzSuperEllipseComputer()
            obj.init();
            obj.solveProblem();
            obj.writeOptimizationVideo();
            obj.printVigdergauzSuperEllipse();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.txi = pi/8;
            obj.rho = 0.8;
        end
        
        function solveProblem(obj)
            s.rho = obj.rho;
            s.txi = obj.txi;
            s.savingFrames = true;
            obj.optimalExponent = OptimalExponentClosestToVigergauz(s);
            obj.optimalExponent.compute();
            obj.xopt   = obj.optimalExponent.qOpt;
            obj.q      = obj.optimalExponent.qOpt;
            obj.error  = obj.optimalExponent.error;    
            obj.frames = obj.optimalExponent.frames;
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
            s.mx = obj.optimalExponent.mx;
            s.my = obj.optimalExponent.my;
            s.txi = obj.txi;
            s.rho = obj.rho;
            s.q   = obj.q;
            stressPrinter = VigdergauzSuperEllipseStressPrinter(s);
            stressPrinter.print();            
        end
        
    end
    
end

