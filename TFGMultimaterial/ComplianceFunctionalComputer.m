classdef ComplianceFunctionalComputer < handle

    properties (Access = public)
        J
        dJ
    end

    properties (Access = private)
        designVariable
        forces
        displ
        energy0
        nMat
        dC
        strain
        tgamma
    end

    methods (Access = public)
        
        function obj = ComplianceFunctionalComputer(cParams)
            obj.init(cParams);
            obj.computeFunction();
            obj.computeGradient();
        end
    end

    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;
            obj.forces = cParams.forces;
            obj.displ = cParams.displacements;
            obj.energy0 = cParams.energy0;
            obj.nMat = cParams.nMat;
            obj.dC = cParams.dC;
            obj.strain = cParams.strain;
            obj.tgamma = cParams.tgamma;
        end

        function computeFunction(obj)
            F = obj.forces;
            U = obj.displ;
        
            obj.J = 0.5*dot(F,U) / obj.energy0; 
        end

        function computeGradient(obj)

            tgamma3 = [obj.tgamma;obj.tgamma;obj.tgamma]; 
            e = obj.strain.*tgamma3;
            for i = 1:obj.nMat
                for j = 1:obj.nMat
                    for z=1:size(obj.dC,3)
                        derTop(z) = e(:,z)'*obj.dC(:,:,z,i,j)*e(:,z);
                    end
                    TD{i,j} = derTop; % lo que abans era TD
                    obj.dJ{i,j} = TD{i,j}/obj.energy0; % adimensionat
                end
            end

        end
    end


end