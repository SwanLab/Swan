classdef TSection < SectionVariablesComputer
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        A
        Iy
        Ix
        dAdtw
        dAdh
        dIydtw
        dIydh
        dIxdtw
        dIxdh
    end
    
    methods (Access = public)
        
        function obj = TSection(cParams)
            obj.init(cParams)
            obj.computeAreaAndInertiaInSymbolic()
        end

       function [Ix,Iy] = computeInertia(obj)
           [h,tw] = obj.getDoubleValue();
           Ix = obj.Ix(h,tw);
           Iy = obj.Iy(h,tw);
        end
        
        function [dIydtw, dIydh, dIxdtw, dIxdh] = computeInertiaDerivative(obj)
           [h,tw] = obj.getDoubleValue();
           dIydtw = obj.dIydtw(h,tw);
           dIydh = obj.dIydh(h,tw);
           dIxdtw = obj.dIxdtw(h,tw);
           dIxdh = obj.dIxdh(h,tw);
        end
        
        function A = computeArea(obj)
           [h,tw] = obj.getDoubleValue();
           A = obj.A(h,tw);

        end

        function dA = computeAreaDerivative(obj)
           [h,tw]=obj.getDoubleValue();
           dAdh = obj.dAdh(h,tw);
           dAdtw = obj.dAdtw(h,tw);
           dA = [dAdh; dAdtw];
        end        
  
    end
    
    methods (Access = private)

        
        function computeAreaAndInertiaInSymbolic(obj)
            tw = sym('tw','real');
            h = sym('h','real');
            b = h;
            tf = tw;
            A  = (b-tw)*tf+h*tw;
            
            Iy = ((h-2*tf)*tw^3/12)+(tf*b^3/12);

            yc = (tw*h^2/2+(b-tw)*tf^2/2);
            Ix1 = tw*h^3/3+(b-tw)*tf^3/3;
            Ix = Ix1-1/A*yc^2;
            
            dAdtw = diff(A,tw);
            dAdh = diff(A,h)
            dIydtw = diff(Iy,tw);
            dIydh = diff(Iy,h);

            dIxdtw = diff(Ix,tw);
            dIxdh = diff(Ix,h);

  
            
            obj.A  = matlabFunction(A);
            obj.dAdtw = matlabFunction(dAdtw);
            obj.dAdh = matlabFunction(dAdh);
            obj.Iy = matlabFunction(Iy);
            obj.Ix = matlabFunction(Ix);
            obj.dIydtw = matlabFunction(dIydtw);
            obj.dIydh  = matlabFunction(dIydh);
            obj.dIxdtw = matlabFunction(dIxdtw);
            obj.dIxdh  = matlabFunction(dIxdh);            
        end

    end
    
end