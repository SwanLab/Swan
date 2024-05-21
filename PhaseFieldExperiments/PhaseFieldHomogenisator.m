classdef PhaseFieldHomogenisator < handle
    
    properties (Access = public)
        E
        nu
        steps
    end   
    
    methods (Access = public)
        
        function obj = PhaseFieldHomogenisator(cParams)
            obj.init(cParams)
        end
        
        function computeIsotropicMaterial(obj,dmg)
            constant = obj.E/(1-obj.nu^2);
            voidLength = linspace(0,1,obj.steps);
            switch dmg
                case 'AT1'
                    phi = voidLength;
                case 'AT2'
                    phi = voidLength.^2;
            end
                      
            mat = zeros(3,3,length(phi));
            mat(1,1,:) = constant*(1-phi).^2 ;
            mat(1,2,:) = constant*obj.nu*(1-phi).^2;
            mat(2,1,:) = constant*obj.nu*(1-phi).^2;
            mat(2,2,:) = constant*(1-phi).^2;
            mat(3,3,:) = constant*(1-obj.nu)*(1-phi).^2 /2;
            
            name = ['IsoMicroDamage',dmg];
            save(name,'mat','phi');
        end
        
        function computeHomogenizedMaterial(obj,voidType,dmgType)
            voidMax = obj.computeLengthMax(voidType);
            phiType = obj.computeDamageModel(voidMax,dmgType);
            voidLength = linspace(0,voidMax,obj.steps);
            mat = zeros(3,3,length(voidLength));
            phi = zeros(length(voidLength),1);
            for i=1:obj.steps
                homogMat = Tutorial02p2FEMElasticityMicro(voidLength(i),voidType,obj.E,obj.nu);
                mat(:,:,i) = homogMat.stateProblem.Chomog;
                
                figure(1)
                cla reset
                homogMat.mesh.plot
            
                phi(i) = phiType(voidLength(i));
            end

            name = [voidType,'MicroDamage',dmgType];
            save(name,'mat','phi');
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.E = cParams.E;
            obj.nu = cParams.nu;
            obj.steps = cParams.steps;
        end
        
        function holeMax = computeLengthMax(~,voidType)
            switch voidType
                case "Circle"
                    holeMax = 0.5;
                case "Square"
                    holeMax = 0.99;
                case "Full"
                    holeMax = 1;
            end
        end
        
        function phiType = computeDamageModel(~,voidMax,dmgType)
            switch dmgType
                case 'Length'
                    phiType = @(x) x/voidMax;
                case 'Volume'
                    phiType = @(x) (x/voidMax)^2;
            end
        end
    end
    
end