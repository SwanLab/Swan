%%  Reads FreeFem++ Data Files into Matlab / Octave
%
%  Author: Chloros2 <chloros2@gmx.de>
%  Created: 2018-12-18
%
% [varargout] = ffreaddata(filename)
%
%  This file is part of the ffmatlib which is hosted at
%  https://github.com/samplemaker/freefem_matlab_octave_plot
%

%% Description
%
%  Reads multidimensional ascii data (complex or real) into the Matlab /
%  Octave workspace. Typically this function is used with the export macros
%  ffExportVh() and ffExportData1() which are located in ffexport.idp.
%

%% Input Parameters
%
%  filename
%

%% Output Parameters
%
%  varargout: double arrays containing the data
%

%% Licence
%
%  Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
%  This program is free software: you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see
%  <https://www.gnu.org/licenses/>.
%

%% Code

function [varargout] = ffreaddata(filename)

    verbose=false;

    if (nargin ~= 1)
        printhelp();
        error('wrong number arguments');
    end
    fid = fopen(filename,'r');
    if fid < 0
        error('cannot open file %s', filename);
    end
    firstLineStr = fgetl(fid);
    frewind(fid);

    %In a FreeFem++ data ascii text file complex numbers may appear as
    %(-5.92176,3.15827) whereas real valued numbers may appear as 0.894621
    %Lexer:
    firstLineCell = strsplit(strtrim(firstLineStr),' ');
    nCols = numel(firstLineCell);
    %Parser:
    formatSpec = '';
    isComplexNo = false(1,nCols);
    for i = 1:nCols
        n1 = numel(strfind(char(firstLineCell(i)),'('));
        n2 = numel(strfind(char(firstLineCell(i)),')'));
        if (n1 > 0)
            %presume complex number
            formatSpec = [formatSpec ' ' '(%f,%f)'];
            isComplexNo(i) = true;
            if (n1 > 1) || (n1 ~= n2)
                error('ffreaddata: parse error complex data');
            end
        else
            %presume real valued number
            formatSpec = [formatSpec ' ' '%f'];
            isComplexNo(i) = false;
        end
    end
    fdata = textscan(fid,strtrim(formatSpec),'Delimiter','\n');
    fclose(fid);

    varargout=cell(1,nCols);
    k = 1;
    for j = 1:nCols
        if isComplexNo(j)
            varargout{j} = fdata{k} + 1i*fdata{k+1};
            k = k + 2;
        else
            varargout{j} = fdata{k};
            k = k + 1;
        end
    end

    if verbose
        fprintf('formatSpec: %s\n',strtrim(formatSpec));
        fprintf('Size of data: (nDof, cols) %ix%i\n',numel(fdata{1}),nCols);
    end

end

function printhelp()
    fprintf('%s\n\n','Invalid call to ffreaddata. Correct usage is:');
    fprintf('%s\n',' -- [varargout] = ffreaddata (filename)');
end
