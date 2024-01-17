classdef HomogenizedMicrostructrueInterpolator < handle

    properties (Access = private)
       Ctensor 
    end
    
    properties (Access = private)
        fileName
    end
    
    methods (Access = public)
        
        function obj = HomogenizedMicrostructrueInterpolator(cParams)
            obj.init(cParams)
            obj.loadConsitutiveTensor();
        end

        function computeConsitutiveTensor(obj,x)
            xR = obj.reshapeDesignVariable(x);
            obj.obtainReferenceConsistutiveTensor(xR);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName = cParams.fileName;
        end
        
        function loadConsitutiveTensor(obj)
            s.fileName = [obj.fileName,'WithAmplificators'];
            v = VademecumVariablesLoader(s);
            v.load();
            obj.Ctensor = v.Ctensor;            
        end

        function obtainReferenceConsistutiveTensor(obj,x)
            mx = x{1};
            my = x{2};
            [c,dc] = obj.Ctensor.compute([mx,my]);
            obj.Cref  = c;
            obj.dCref = permute(dc,[1 2 4 3 5]);
        end        

        function m = createMaterial(obj,C)
            s.type    = 'ANISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.ndim;
            s.constiutiveTensor = C;            
            s.bulk    = kappa;
            m = Material.create(s);   
        end
        
    end
    
    methods (Access = private, Static)
        
        function xR = reshapeDesignVariable(x)
            xV = x.value;
            nV = x.nVariables;
            nx = length(xV)/nV;
            xR = cell(nV,1);
            for ivar = 1:nV
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xR{ivar} = xV(i0:iF);
            end
        end
        
    end    
    
end