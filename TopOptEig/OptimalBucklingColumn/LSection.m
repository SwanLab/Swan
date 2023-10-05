classdef LSection < SectionVariablesComputer
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        A
        Iy
        Ix
        dIydtw
        dIydh
        dIxdtw
        dIxdh
    end
    
    methods (Access = public)
        
        function obj = LSection(cParams)
            obj.init(cParams)
            obj.computeAreaAndInertiaInSymbolic()
        end

       function [Ix,Iy] = computeInertia(obj)
           [h,tw]=obj.getDoubleValue();
           Ix = obj.Ix(h,tw);
           Iy = obj.Iy(h,tw);
        end
        
        function dI = computeInertiaDerivative(obj)
           [h,tw]=obj.getDoubleValue();

        end
        
        function A = computeArea(obj)
           [h,tw]=obj.getDoubleValue();

        end

        function dA = computeAreaDerivative(obj)
           [h,tw]=obj.getDoubleValue();

        end        
  
    end
    
    methods (Access = private)

        
        function computeAreaAndInertiaInSymbolic(obj)
            tw = sym('tw','real');
            h = sym('h','real');
            b = h;
            tf = tw;
            A  = tw*h + tf*b - tf*tw;
            xC = 1/(2*A)*(h*tw^2 + b^2*tf-tw^2*tf);
            Iy0 = 1/3*(h*tw^3+b^3*tf-tw^3*tf);
            Iy = Iy0 - A*xC^2;

            dIydtw = diff(Iy,tw);
            dIydh = diff(Iy,h);


            yC = 1/(2*A)*(h^2*tw + b*tf^2-tw*tf^2);
            Ix0 = 1/3*(h^3*tw+b*tf^3-tw*tf^3);
            Ix = Ix0 - A*yC^2;

            dIxdtw = diff(Ix,tw);
            dIxdh = diff(Ix,h);

  
            
            obj.A  = matlabFunction(A);
            obj.Iy = matlabFunction(Iy);
            obj.Ix = matlabFunction(Ix);
            obj.dIydtw = matlabFunction(dIydtw);
            obj.dIydh  = matlabFunction(dIydh);
            obj.dIxdtw = matlabFunction(dIxdtw);
            obj.dIxdh  = matlabFunction(dIxdh);            
        end

    end
    
end