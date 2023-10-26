classdef LSection < SectionVariablesComputer
      
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
        ind
    end
    
    methods (Access = public)
        
        function obj = LSection(cParams)
            obj.init(cParams)
            obj.computeAreaAndInertiaInSymbolic()
        end

       function I= computeInertia(obj)
           [h,tw] = obj.getDoubleValue();
           Ix = obj.Ix(h,tw);
           Iy = obj.Iy(h,tw);
           [I,index] = min([Ix,Iy],[],2);
           obj.ind = index;
        end
        
        function dI = computeInertiaDerivative(obj)
           [h,tw] = obj.getDoubleValue();
           dIydtw = obj.dIydtw(h,tw);
           dIydh = obj.dIydh(h,tw);
           dIxdtw = obj.dIxdtw(h,tw);
           dIxdh = obj.dIxdh(h,tw);
           
           dIdtwT = [dIxdtw,dIydtw];
           dIdhT = [dIxdh,dIydh];

           idx = sub2ind(size(dIdtwT), [1:size(dIdtwT,1)]', obj.ind);
           dI = [dIdtwT(idx); dIdhT(idx)];


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
            A  = tw*h + tf*b - tf*tw;
            xC = 1/(2*A)*(h*tw^2 + b^2*tf-tw^2*tf);
            Iy0 = 1/3*(h*tw^3+b^3*tf-tw^3*tf);
            Iy = Iy0 - A*xC^2;
            
            dAdtw = diff(A,tw);
            dAdh = diff(A,h);
            dIydtw = diff(Iy,tw);
            dIydh = diff(Iy,h);


            yC = 1/(2*A)*(h^2*tw + b*tf^2-tw*tf^2);
            Ix0 = 1/3*(h^3*tw+b*tf^3-tw*tf^3);
            Ix = Ix0 - A*yC^2;

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