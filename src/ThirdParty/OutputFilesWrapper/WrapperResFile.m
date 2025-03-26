classdef WrapperResFile < FileReader
    
    properties (Access = private)
       dataBase
       stress
       exponent
       alpha
       lineNumber
       linesEndValue
    end
    
    properties (Access = private)
       nElem
       nNodes
       dimension       
    end
    
    methods (Access = public)
        
        function obj = WrapperResFile(cParams)
            obj.init(cParams);            
        end
        
        function read(obj)
            obj.openFile();
            obj.readStress();
            obj.readStrain();  
            obj.readDisplacement();  
            obj.readStress();
            obj.readStrain();  
            obj.readDisplacement();              
            obj.readAmplifiedStressNormP();
            %obj.readCompliance();
            obj.readAngle();  
            obj.readAngle(); 
            obj.readAmplifiedStressNormP();            
            obj.readDensity();  
            obj.readScalarInNodes();
            obj.readScalarInNodes();
            obj.readScalarInGauss();            
            obj.readScalarInNodes();
            obj.readScalarInNodes();
            obj.closeFile();
        end
        
        function d = getDataBase(obj)
            d = obj.dataBase;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePath  = cParams.filePath;
            obj.nElem     = cParams.nElem;
            obj.nNodes    = cParams.nNodes;
            obj.dimension = cParams.dimension;
            obj.lineNumber = 14;    
            obj.linesEndValue = 2;            
        end
        
        function readStress(obj)
            obj.readTensorFieldInGauss()
        end
        
        function readStrain(obj)
            obj.readTensorFieldInGauss()
        end
        
        function readDisplacement(obj)
            obj.readVectorFieldInNodes();
        end
        
        function readCompliance(obj)
            obj.readScalarInGauss();            
        end
        
        function readAmplifiedStressNormP(obj)
            obj.readScalarInGauss();            
        end
                
        function readAngle(obj)
            obj.readVectorFieldInGauss();              
        end
        
        function readDensity(obj)
            obj.readScalarInGauss();            
        end
        
        function readLevelSet(obj)
            obj.readScalarInNodes();
        end
               
        
        function readVectorFieldInGauss(obj)
            obj.readVectorInGauss(obj.dimension)
        end
        
        function readVectorFieldInNodes(obj)
            obj.readVectorInNodes(obj.dimension)
        end        
        
        function readTensorFieldInGauss(obj)
            switch obj.dimension                
                case 2
                    nComp = 3;
                case 3
                    nComp = 6;
            end
            obj.readVectorInGauss(nComp);            
        end
        
        function readScalarInGauss(obj)
            nLines = obj.nElem;
            linesJump = 2;
            nComp = 1;
            obj.readField(obj.lineNumber,nLines,linesJump,nComp);
            obj.updateLineNumber(nLines + linesJump + obj.linesEndValue);
        end
        
        function readVectorInGauss(obj,nComp)
            nLines = obj.nElem;
            linesJump = 3;
            obj.readField(obj.lineNumber,nLines,linesJump,nComp);
            obj.updateLineNumber(nLines + linesJump + obj.linesEndValue);           
        end        
        
        function readVectorInNodes(obj,nComp)
            nLines = obj.nNodes;
            linesJump = 3;
            obj.readField(obj.lineNumber,nLines,linesJump,nComp);
            obj.updateLineNumber(nLines + linesJump + obj.linesEndValue);           
        end             
        
        function readScalarInNodes(obj)
            nLines = obj.nNodes;
            linesJump = 2;
            nComp = 1;
            obj.readField(obj.lineNumber,nLines,linesJump,nComp);
            obj.updateLineNumber(nLines + linesJump + obj.linesEndValue);            
        end
        
        function [fName,fValue] = readField(obj,lineNumber,nLines,linesJump,nComp)
            fName = obj.readFieldName(lineNumber);
            lineNumber = lineNumber + linesJump;            
            fValue = obj.readlFieldValues(lineNumber,nLines,nComp);
            obj.dataBase.(fName) = fValue;                        
        end        
        
        function name = readFieldName(obj,lineNumber)
            lineText = obj.readLines(lineNumber,1); 
            lSplit = split(lineText{1});
            name = lSplit(2);
            name = obj.eraseQuationMarks(name{1});
        end        
        
        function fValue = readlFieldValues(obj,lineNumber,nLines,nComp)
            lines = obj.readLines(lineNumber,nLines);
            namesSplit = split(lines{:});
            lines = str2double(namesSplit);
            index = (1:nComp) + 1;
            fValue = lines(:,index); 
        end
        
        function tLine = readLines(obj,lineNumber,nLines)
            fseek(obj.fid,0,'bof');
            tLine = textscan(obj.fid,'%s',nLines,'delimiter','\n', 'headerlines',lineNumber-1);
        end
        
        function updateLineNumber(obj,lines)
            obj.lineNumber = obj.lineNumber + lines;
        end
        
    end
    
    methods (Access = private, Static)
        
        function str = eraseQuationMarks(str)
            str = erase(str,"""");
        end
        
    end
    
end