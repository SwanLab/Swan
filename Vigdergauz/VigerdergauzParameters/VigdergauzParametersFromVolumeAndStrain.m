classdef VigdergauzParametersFromVolumeAndStrain < handle
    
    properties (Access = public)
        parameters
    end
    
    properties (Access = private)
        ax
        ay
        cx
        cy
        strain
        volume
        mu1
        mu0
        k1
        k0
    end
    
    methods (Access = public)
        
        function obj = VigdergauzParametersFromVolumeAndStrain(cParams)
            obj.init(cParams);
        end
        
        function compute(obj)
            obj.computeAxAy();
            obj.computeVigergauzParameters();            
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.strain = cParams.strain;
            obj.volume = cParams.volumeMicro;
            obj.cx     = cParams.cx;
            obj.cy     = cParams.cy;            
            k = @(E,nu) E/(2*(1-nu));
            mu = @(E,nu) E/(2*(1+nu));
            obj.mu0 = mu(cParams.E0,cParams.nu0);
            obj.k0  = k(cParams.E0,cParams.nu0);
            obj.mu1 = mu(cParams.E1,cParams.nu1);            
            obj.k1  = k(cParams.E1,cParams.nu1);
        end
        
        function computeAxAy(obj)
            s.strain = obj.strain;
            s.volume = obj.volume;
            s.cx     = obj.cx;
            s.cy     = obj.cy;
            s.mu0    = obj.mu0;
            s.mu1    = obj.mu1;
            s.k0    = obj.k0;
            s.k1    = obj.k1;            
            axay = AxAyComputerFromVolumeAndStrain(s);
            [obj.ax,obj.ay] = axay.compute();                        
        end
        
        function computeVigergauzParameters(obj)
            s.ax = obj.ax;
            s.ay = obj.ay;
            s.cx = obj.cx;
            s.cy = obj.cy;
            s.type = 'AxAndAy';
            obj.parameters = VigdergauzParameters.create(s);
        end         
        
    end
    
end