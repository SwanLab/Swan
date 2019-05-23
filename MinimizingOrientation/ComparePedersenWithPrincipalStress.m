classdef ComparePedersenWithPrincipalStress < handle
    
    properties (Access = private)
        Ctensor
        strain
        principalDirection
        optimalAngle
        minimalCompOrientation
    end
    
    methods (Access = public)
        
        function obj = ComparePedersenWithPrincipalStress()
            obj.init();
            obj.computePrincipalStress();
            obj.computeMinimalComplianceOrientation();
            obj.compareResults();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.createStrain();
            %obj.createOrthotropicConstitutiveTensor();
            obj.createParticularOrthotropicConstitutiveTensor
            %obj.createIsotropicConstitutiveTensor
        end
        
        function createStrain(obj)
            obj.strain = [1;0.3;0];
            %obj.strain = [2.8048843;-0.40821;0];
            %obj.strain = rand(3,1);
        end
        
        function createOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = rand(1,1);
            C(2,2) = rand(1,1);
            C(1,2) = obj.randBetweenBounds(-0.1,0.1);
            C(2,1) = C(1,2);
            C(3,3) = obj.randBetweenBounds(0,1);
            obj.Ctensor = C;
            det(C)
        end
        
        function createParticularOrthotropicConstitutiveTensor(obj)
            C = zeros(3,3);
            C(1,1) = 0.374827;
            C(2,2) = 0.864285;
            C(1,2) = 0.1257865;
            C(2,1) = C(1,2);
            C(3,3) = 0.13265;
            obj.Ctensor = C;
            det(C)
        end        
        
        function createIsotropicConstitutiveTensor(obj)
            E = obj.randBetweenBounds(0.01,1);
            nu = obj.randBetweenBounds(0,0.49);
            C = zeros(3,3);
            C(1,1) = 1;
            C(2,2) = 1;
            C(1,2) = nu;
            C(2,1) = C(1,2);
            C(3,3) = (1-nu)/2;
            C = E*C;
            obj.Ctensor = C;
        end
        
        function computePrincipalStress(obj)
            cParams.strain  = obj.strain;
            cParams.Ctensor = obj.Ctensor;
            ps = PrincipalStressDirections(cParams);
            ps.compute();
            obj.principalDirection = ps;
        end
        
        function computeMinimalComplianceOrientation(obj)
            p = Pedersen();
            p.compute(obj.strain,obj.Ctensor);
            obj.minimalCompOrientation = p;
            obj.optimalAngle = p.optimalAngle;            
        end
        
        function compareResults(obj)
            a = obj.optimalAngle;
            Rp = [cos(a) -sin(a); sin(a) cos(a)];
            Rd = obj.principalDirection.principalStressDir;
            stress = obj.principalDirection.stress;
            Rsig = obj.minimalCompOrientation.RotationMatrix;
            norm(Rp - Rd)
            ad = acos(Rd(1,1));
            comp = obj.minimalCompOrientation.compliance;
            xp = linspace(0,2*pi,100);
            plot(xp*180/pi,comp(xp))            
            hold on
            plot(a*180/pi,comp(a),'b+');
            plot(ad*180/pi,comp(ad),'r+');
            plot((ad+pi/2)*180/pi,comp(ad+pi/2),'g+');
        end
        
    end
    
    methods (Access = private, Static)
        
        function c = randBetweenBounds(a,b)
            c = (b-a).*rand(1,1) + a;
        end
        
    end
    
    
    
end