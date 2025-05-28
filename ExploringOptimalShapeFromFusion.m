classdef ExploringOptimalShapeFromFusion < handle

    properties (Access = private)
        filename
        mesh       
        young
        poisson   
        material
        physicalProblem
        compliance
        volume
    end

    methods (Access = public)

        function obj = ExploringOptimalShapeFromFusion()
            obj.init()
            obj.createMesh();
            obj.createVolume();
        end

    end

    methods (Access = private)

        function init(obj)
            close all;
        end

        function createMesh(obj)
            file = 'BEAM_3D_FRAME';
            obj.filename = file;
            a.fileName = file;
            s = FemDataContainer(a);
            obj.mesh = s.mesh;
        end


        function createVolume(obj)
            volume = obj.mesh.computeVolume();
            InitialVolume = 10;
            fractionVolume = volume/InitialVolume;
        end        


        function computeElasticProperties(obj)
            E  = 1;
            nu = 1/3;
            obj.young   = ConstantFunction.create(E,obj.mesh);
            obj.poisson = ConstantFunction.create(nu,obj.mesh);
        end

        function createMaterial(obj)
            s.type    = 'ISOTROPIC';
            s.ptype   = 'ELASTIC';
            s.ndim    = obj.mesh.ndim;
            s.young   = obj.young;
            s.poisson = obj.poisson;
            tensor    = Material.create(s);
            obj.material = tensor;
        end


      

    end


end