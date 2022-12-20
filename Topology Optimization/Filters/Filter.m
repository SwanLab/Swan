classdef Filter < handle

    properties (Access = protected)
        field
        fieldM
    end

    properties (Access = protected)
        mesh
        quadratureOrder
        femSettings
        LHStype
    end

    methods(Access = public, Static)

        function obj = create(cParams)
            f = FilterFactory();
            obj = f.create(cParams);
        end

    end

    methods(Access = protected)

        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.createField();
            obj.createFieldMass();
            obj.femSettings = cParams.femSettings;
            if isfield(cParams.femSettings,'LHStype')
                obj.LHStype = cParams.femSettings.LHStype;
            else
                obj.LHStype = 'DiffReactNeumann';
            end
        end

        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'LINEAR';
            obj.field = Field(s);
        end

        function createFieldMass(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.fieldM = Field(s);
        end

    end

end