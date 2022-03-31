classdef EigenValueAndVectorComputerPrecomputed < EigenValueAndVectorComputer
    
    methods (Access = public)
        
        function obj = EigenValueAndVectorComputerPrecomputed(cParams)
            obj.init(cParams);
            obj.computeEigenVectorAndValueFunction();
        end
        
    end
    
    methods (Access = private)
        
        function computeEigenVectorAndValueFunction(obj)
            switch obj.ndim
                case 2
                    obj.computeEigenValueFunction2D();
                    obj.computeEigenVectorFunction2D();
                case 3
                    obj.computeEigenValueFunction3D();
                    obj.computeEigenVectorFunction3D();
            end
        end
        
        function computeEigenVectorFunction2D(obj)
            eV{1,1} = @(S1_1,S2_1,S2_2)-(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).*1.0./sqrt(abs(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).^2+1.0);
            eV{1,2} = @(S1_1,S2_1,S2_2)-1.0./sqrt(abs(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).^2+1.0);
            eV{2,1} = @(S1_1,S2_1,S2_2)1.0./sqrt(abs(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).^2+1.0);
            eV{2,2} = @(S1_1,S2_1,S2_2)-(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).*1.0./sqrt(abs(S2_2./S2_1-(S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0)./S2_1).^2+1.0);
            obj.eigenVectorFunction = eV;
        end
       
        function computeEigenValueFunction2D(obj)
            eV{1} = @(S1_1,S2_1,S2_2)S1_1./2.0+S2_2./2.0-sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0;
            eV{2} = @(S1_1,S2_1,S2_2)S1_1./2.0+S2_2./2.0+sqrt(S1_1.*S2_2.*-2.0+S1_1.^2+S2_1.^2.*4.0+S2_2.^2)./2.0;
            obj.eigenValueFunction = eV;
        end
        
        
        function computeEigenVectorFunction3D(obj)
            d = load('EigenValueAndVectorComputer3D.mat');
            eV = d.eigenComputer;
            obj.eigenVectorFunction = eV.eigenVectorFunction;
        end
       
        function computeEigenValueFunction3D(obj)
            d = load('EigenValueAndVectorComputer3D.mat');
            eV = d.eigenComputer;
            obj.eigenValueFunction = eV.eigenValueFunction;
        end
        
    end
    
end