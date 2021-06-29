classdef FreeFemInputPrinter < FreeFemFilePrinter
    
    
    properties (Access = private)
        inputLines
        beforeInputLines
        afterInputLines
        
        startingLines
        endingLines
        inputData
    end
    
    methods (Access = public)
        
        function setInputData(obj,d)
            obj.inputData = d;
        end
    end
    
    methods (Access = protected)
        
        function printingLines(obj)
            obj.obtainBeforeAfterAndInputLines();
            obj.printBeforeInputLines();
            obj.printInputLines();
            obj.printAfterInputLines();
        end
        
    end
    
    methods (Access = private)
        
        function obtainBeforeAfterAndInputLines(obj)
            nlines = numel(obj.linesToPrint);
            [sL,eL] = obj.obtainStartingEndingInputLines();
            obj.beforeInputLines = (1:sL-1);
            obj.inputLines = (sL:eL);
            obj.afterInputLines = (eL+1:nlines);
        end
        
        function [sL,eL] = obtainStartingEndingInputLines(obj)
            nlines = numel(obj.linesToPrint);
            obj.startingLines = false(nlines,1);
            obj.endingLines   = false(nlines,1);
            for iline = 1:nlines
                line = obj.linesToPrint{iline};
                obj.startingLines(iline) = obj.isInputStarting(line);
                obj.endingLines(iline)   = obj.isInputEnding(line);
            end
            sL = obj.obtainStartingLine();
            eL = obj.obtainEndingLine();
        end
        
        function sL = obtainStartingLine(obj)
            sLines = obj.startingLines;
            sLines = obj.booleanIndex2integerIndex(sLines);
            firstSL = sLines(1);
            sL = firstSL + 1;
        end
        
        function eL = obtainEndingLine(obj)
            eLines = obj.endingLines;
            eLines = obj.booleanIndex2integerIndex(eLines);
            firstEL = eLines(1);
            eL = firstEL(1) - 1;
        end
        
        function printBeforeInputLines(obj)
            lines = obj.beforeInputLines;
            obj.printNonInputLines(lines)
        end
        
        function printAfterInputLines(obj)
            lines = obj.afterInputLines;
            obj.printNonInputLines(lines)
        end
        
        function printNonInputLines(obj,nonInputLines)
            nLines = length(nonInputLines);
            for iline = 1:nLines
                lineNumber = nonInputLines(iline);
                line = obj.linesToPrint{lineNumber};
                obj.printLine(line);
            end
        end
        
        function printInputLines(obj)
            obj.printRealInputData();
            obj.printStringInputData();
        end
        
        function printRealInputData(obj)
            inVar = obj.inputData.reals;
            for iData = 1:numel(inVar)
                varName  = inVar{iData}{1,1};
                varValue = num2str(inVar{iData}{1,2});
                line = ['real ',varName,' = ',varValue,';'];
                obj.printLine(line);
            end   
        end
        
        function printStringInputData(obj)
            inVar = obj.inputData.strings;
            for iData = 1:numel(inVar)
                varName  = inVar{iData}{1,1};
                varValue = ['"',inVar{iData}{1,2},'"'];
                line = ['string ',varName,' = ',varValue,';'];
                obj.printLine(line);
            end               
        end
        
        function printInputData(obj)
            inVar = obj.inputData.reals;
            for iData = 1:numel(inVar)
                varName  = inVar{iData}{1,1};
                varValue = num2str(inVar{iData}{1,2});
                line = ['real ',varName,' = ',varValue,';'];
                obj.printLine(line);
            end            
        end
                
        function printInputLine(obj)
            fprintf(obj.fileID,'%s\r\n',obj.line);
        end
        
        function printLine(obj,line)
            fprintf(obj.fileID,'%s\r\n',line);
        end
        
    end
    
    methods (Access = private, Static)
        
        function itIs = isInputStarting(line)
            itIs = strcmp(line,'// --- Input ---');
        end
        
        function itIs = isInputEnding(line)
            itIs = strcmp(line,'// --- End Input ---');
        end
        
        function index = booleanIndex2integerIndex(index)
            allValues = 1:length(index);
            index = allValues(index);
        end
        
    end
    
end