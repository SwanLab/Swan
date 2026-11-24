classdef SidesLatticeDatasetGenerator < handle

    properties (Access = private)
        generator
        transformation
        meshSettings
        motifSettings
        hole
        homogenizer
    end

    methods (Access = public)

        function obj = SidesLatticeDatasetGenerator(cParams)
            obj.generator      = LatticeDatasetGenerator(cParams);
            % obj.transformation = HoleSidesTransformation();
            obj.transformation = SidesLatticeDatasetGenerator.createTransformation(cParams);
            obj.meshSettings   = cParams.mesh;
            obj.motifSettings  = cParams.motif;
            obj.hole           = cParams.hole;
            obj.homogenizer    = MicroHomogenizer(cParams.material);
        end

        function [C,rho] = computeHomogenizationAt(obj,a,sides)
            T            = obj.transformation.compute(a);
            periodicCell = obj.createCell(T);
            mesh         = periodicCell.mesh;
            nSides       = size(sides,1);
            C            = zeros([repmat(mesh.ndim,1,4),nSides]);
            rho          = zeros(nSides,1);
            for i = 1:nSides
                m1 = sides(i,1);
                m2 = sides(i,2);
                if m1*m2 == 0
                    % [C(:,:,:,:,i),rho(i)] = obj.generator.computeHomogenizationAtScale(0,[a,1]);
                    [C(:,:,:,:,i),rho(i)] = obj.generator.computeHomogenizationAtScale(0,[0,1]);
                else
                    motif  = obj.createMotif(m1,m2);
                    dens   = motif.computeDensity(mesh,T,sqrt(m1*m2));
                    rho(i) = motif.computeVolumeFraction(mesh,dens);
                    C(:,:,:,:,i) = obj.homogenizer.compute(periodicCell,dens);
                end
            end
        end

        function div = computeDivisions(obj,a)
            T       = obj.transformation.compute(a);
            [~,div] = PeriodicCell.discretize(T,obj.meshSettings.divUnit,obj.meshSettings.strategy);
        end

        function [X,Y,isSolid,T] = sampleMicrostructure(obj,a,m1,m2,nCells,resolution)
            T       = obj.transformation.compute(a);
            s       = linspace(0,nCells,resolution*nCells+1);
            [S1,S2] = ndgrid(s,s);
            x       = T*[S1(:)';S2(:)'];
            if m1*m2 == 0
                phi = ones(1,size(x,2));
            else
                motif = obj.createMotif(m1,m2);
                phi   = motif.evaluateLevelSet(x,T,sqrt(m1*m2));
            end
            X       = reshape(x(1,:),size(S1));
            Y       = reshape(x(2,:),size(S1));
            isSolid = reshape(phi > 0,size(S1));
        end

        function sym = getSymmetry(obj,name)
            sym = obj.transformation.getSymmetry(name);
        end

        function y = tensorToVector(obj,C)
            y = obj.generator.tensorToVector(C);
        end

        function names = getComponentNames(obj)
            names = obj.generator.getComponentNames();
        end

    end

    methods (Static, Access = public)

        function rho = computeRho(m1,m2)
            rho = 1 - m1.*m2;
        end

        function xi = computeXi(m1,m2)
            xi = atan2(m1,m2);
        end

        function transformation = createTransformation(cParams)
            if isfield(cParams,'cellTransformation')
                transformation = createLatticeTransformation(cParams);
            else
                transformation = HoleSidesTransformation();
            end
        end

    end

    methods (Access = private)

        function c = createCell(obj,T)
            s   = obj.meshSettings;
            s.T = T;
            c   = PeriodicCell(s);
        end

        function motif = createMotif(obj,m1,m2)
            s             = obj.motifSettings;
            s.type        = 'Rectangle';
            s.aspectRatio = sqrt(m1/m2);
            motif         = CellMotif(s);
        end

    end

end