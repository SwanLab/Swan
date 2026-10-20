classdef PlotterMultiField < handle

    properties (Access = private)
        patchHandles
        meshes
        designVariable
        figures
    end

    methods (Access = public)

        function obj = PlotterMultiField(cParams)
            obj.init(cParams);
            obj.createFigures();
        end

        function plot(obj)
            fields = obj.designVariable.fields;
            nFields = numel(fields);

            for iField = 1:nFields
                fieldI = fields{iField};
                funp0 = fieldI.project('P0');
                fieldElem = squeeze(funp0.fValues);

                set(...
                    obj.patchHandles{iField},...
                    'FaceVertexAlphaData',fieldElem,...
                    'FaceAlpha','flat');
            end
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.designVariable = cParams.designVariable;

            fields = obj.designVariable.fields;
            nFields = numel(fields);

            obj.meshes = cell(nFields,1);

            for iField = 1:nFields
                obj.meshes{iField} = fields{iField}.mesh;
            end

            obj.patchHandles = cell(nFields,1);
            obj.figures = gobjects(nFields,1);
        end

        function createFigures(obj)
            nFields = numel(obj.meshes);

            for iField = 1:nFields
                mesh = obj.meshes{iField};

                obj.figures(iField) = figure;
                set(gcf,'Pointer','arrow','NumberTitle','off');

                hold on
                axis off
                axis equal

                ax = gcf().Children;

                obj.patchHandles{iField} = patch(...
                    ax,...
                    'Faces',mesh.connec,...
                    'Vertices',mesh.coord,...
                    'EdgeColor','none',...
                    'LineStyle','none',...
                    'FaceLighting','none',...
                    'AmbientStrength',.75);

                title(['Field ',num2str(iField)]);
            end
        end

    end

end