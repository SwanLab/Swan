classdef SimulationPrinter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        filename
        outputFile
        stepsFiles
    end

    properties (Access = private)
    end
    
    methods (Access = public)
        
        function obj = SimulationPrinter(cParams)
            obj.init(cParams);
            obj.openFile();
        end

        function appendStep(obj, filename)
            obj.stepsFiles{end+1} = filename;
        end

        function print(obj)
            obj.createSimulation();
            obj.saveFile();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj, cParams)
            obj.filename = cParams.filename;
        end

        function openFile(obj)
            fullfile = strcat(obj.filename, '.pvd');
            obj.outputFile = fopen(fullfile,'w');
        end
        
        function createSimulation(obj)
            docNode = com.mathworks.xml.XMLUtils.createDocument('VTKFile');
            % Creating VTKFile
            fileN = docNode.getDocumentElement;
            fileN.setAttribute('type', 'Collection');
            fileN.setAttribute('version', '0.1');
            fileN.setAttribute('byte_order', 'LittleEndian');

            % Creating Collection
            collN = docNode.createElement('Collection');
            fileN.appendChild(collN);
            
            % Add time steps
            for iFun = 1:numel(obj.stepsFiles)
                n = obj.createDataSetNode(docNode, iFun);
                collN.appendChild(n);
            end

            text = xmlwrite(docNode);
            fprintf(obj.outputFile, text);
        end
        
        function saveFile(obj)
        end

        function n = createDataSetNode(obj, docNode, iFun)
            n = docNode.createElement('DataSet');
            n.setAttribute('timestep', num2str(iFun));
            n.setAttribute('group', '');
            n.setAttribute('part', num2str(0));
            n.setAttribute('file', [obj.stepsFiles{iFun}, '.vtu']);
            n.setAttribute('name', obj.filename);
        end
        
    end
    
end

