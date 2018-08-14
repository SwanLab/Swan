function varargout = imclipboard(clipmode, varargin)
%IMCLIPBOARD Copy and paste images to and from system clipboard.
%
%   IMCLIPBOARD('copy', IMDATA) sets the clipboard content to the image
%   represented by IMDATA. IMDATA must be MxN grayscale (double, uint8,
%   uint16), MxN black and white (logical), MxNx3 true color (double,
%   uint8, uint16)
%
%   IMCLIPBOARD('copy', X, MAP) sets the clipboard content to the image
%   data represented by indexed image X with colormap MAP. X must be MxN
%   matrix (double, uint8, uint16) and MAP must be Px3 (double).
%
%   IMDATA = IMCLIPBOARD('paste') returns the current image content in the
%   clipboard as a true color image (MxNx3 uint8).
%
%   [X, MAP] = IMCLIPBOARD('paste') returns the current image content in
%   the clipboard as an indexed color image.
%
%   IMCLIPBOARD('paste') displays the image in a new figure window.
%
%   [...] = IMCLIPBOARD('paste', FILENAME) saves the image as FILENAME. FILENAME
%   must be a name to one of the following image formats: JPG, GIF, BMP,
%   PNG, TIF.
%
%   Note: IMCLIPBOARD requires Java on all platforms.
%
%   Example:
%       im = imread('peppers.png');
%       imclipboard('copy', im);
%       % paste into a paint program
%
%       im2 = imclipboard('paste');
%
%   See also CLIPBOARD.

% Jiro Doke
% Copyright 2010 The MathWorks, Inc.
% Sept 12, 2010

narginchk(1, 3);
error(javachk('awt', 'IMCLIPBOARD'));

% Import necessary Java classes
import java.awt.Toolkit
import java.awt.image.BufferedImage
import java.awt.datatransfer.DataFlavor

% Get System Clipboard object (java.awt.Toolkit)
cb = Toolkit.getDefaultToolkit.getSystemClipboard();

clipmode = validatestring(clipmode, {'copy', 'paste'}, mfilename, 'CLIPMODE');

switch clipmode
    case 'copy'
        nargoutchk(0, 0);

        % Add java class (ImageSelection) to the path
        if ~exist('ImageSelection', 'class')
            javaaddpath(fileparts(which(mfilename)), '-end');
        end

        data = validateCopyInput(varargin{:});
        
        % Get image size
        ht = size(data, 1); wd = size(data, 2);
        
        % Convert to Blue-Green-Red format
        data = data(:, :, [3 2 1]);
        
        % Convert to 3xWxH format
        data = permute(data, [3, 2, 1]);
        
        % Append Alpha data (not used)
        data = cat(1, data, 255*ones(1, wd, ht, 'uint8'));
        
        % Create image buffer
        imBuffer = BufferedImage(wd, ht, BufferedImage.TYPE_INT_RGB);
        imBuffer.setRGB(0, 0, wd, ht, typecast(data(:), 'int32'), 0, wd);
        
        % Create ImageSelection object
        %    % custom java class
        imSelection = ImageSelection(imBuffer);
        
        % Set clipboard content to the image
        cb.setContents(imSelection, []);
                
    case 'paste'
        nargoutchk(0, 2);

        filename = validatePasteInput(varargin{:});
        
        try
            % Attempt to retrieve image data from system clipboard. If there is no
            % image data, it will throw an exception.
            imBuffer = cb.getData(DataFlavor.imageFlavor);
        catch %#ok<CTCH>
            disp('No image data in clipboard');
            imBuffer = [];
        end
        
        if isempty(imBuffer)
            im = [];
        else
            im = imBuffer.getRGB(0, 0, imBuffer.getWidth, imBuffer.getHeight, [], 0, imBuffer.getWidth);
            % "im" is an INT32 array, where each value contains 4 bytes of information
            % (Blue, Green, Red, Alpha). Alpha is not used.
            
            % type cast INT32 to UINT8 (--> 4 times as many elements)
            im = typecast(im, 'uint8');
            
            % Reshape to 4xWxH
            im = reshape(im, [4, imBuffer.getWidth, imBuffer.getHeight]);
            
            % Remove Alpha information (4th row) because it is not used
            im = im(1:3, :, :);
            
            % Convert to HxWx3 array
            im = permute(im, [3 2 1]);
            
            % Convert color space order to R-G-B
            im = im(:, :, [3 2 1]);
        end
        
        if nargout == 1
            varargout{1} = im;
        elseif nargout == 2
            if isempty(im)
                varargout{1} = [];
                varargout{2} = [];
            else
                % Convert to indexed image (getting full possible color space)
                [varargout{1}, varargout{2}] = rgb2ind(im, 65536);
            end
        else
            if nargin == 1 && ~isempty(im)
                figure;
                imshow(im);
            end
        end
        if nargin == 2 && ~isempty(im)
            [~, ~, e] = fileparts(filename);
            if strcmpi(e, '.gif')
                % GIF files require indexed image
                [x, map] = rgb2ind(im, 256);
                imwrite(x, map, filename);
            else
                imwrite(im, filename);
            end
        end
        
end

end

%% Helper Function
function data = validateCopyInput(varargin)

if nargin == 0
    error('imclipboard:NotEnoughArguments', ...
        'For ''copy'', the second argument must be an image array');
end

if nargin == 2  % X and colormap
    % X
    validateattributes(varargin{1}, {'double', 'uint8', 'uint16'}, {'2d', 'real', 'nonnegative'}, mfilename, 'X');
    if isa(varargin{1}, 'double')
        validateattributes(varargin{1}, {'double'}, {'>=', 0, '<=', 1}, mfilename, 'X');
    end
    
    % MAP
    validateattributes(varargin{2}, {'double'}, {'size', [nan 3], '>=', 0, '<=', 1}, mfilename, 'MAP');
    data = ind2rgb(varargin{1}, varargin{2});
    data = uint8(255*data);

else % nargin == 1 % Grayscale, True Color, Black and White
        validateattributes(varargin{1}, ...
        {'double', 'logical', 'uint8', 'uint16'}, ...
        {'real', 'nonnegative'}, mfilename, 'IMDATA');

    switch class(varargin{1})
        case 'double'
            validateattributes(varargin{1}, {'double'}, {'>=', 0, '<=', 1}, mfilename, 'IMDATA');
            if ndims(varargin{1}) == 2  % grayscale
                data = uint8(255*varargin{1});
                data = cat(3, data, data, data);
            elseif ndims(varargin{1}) == 3
                validateattributes(varargin{1}, {'double'}, {'size', [nan nan 3]}, mfilename, 'IMDATA');
                data = uint8(255*varargin{1});
            else
                error('imclipboard:InvalideDoubleImage', ...
                    'Double image data must be 2 or 3 dimensions');
            end
            
        case 'logical'
            validateattributes(varargin{1}, {'logical'}, {'2d'}, mfilename, 'IMDATA');
            data = uint8(255*varargin{1});
            data = cat(3, data, data, data);
            
        case 'uint8'
            if ndims(varargin{1}) == 2
                data = cat(3, varargin{1}, varargin{1}, varargin{1});
            elseif ndims(varargin{1}) == 3
                validateattributes(varargin{1}, {'uint8'}, {'size', [nan nan 3]}, mfilename, 'IMDATA');
                data = varargin{1};
            else
                error('imclipboard:InvalideUINT8Image', ...
                    'UINT8 image data must be 2 or 3 dimensions');
            end
            
        case 'uint16'
            if ndims(varargin{1}) == 2
                data = uint8(varargin{1}*(double(intmax('uint8'))/double(intmax('uint16'))));
                data = cat(3, data, data, data);
            elseif ndims(varargin{1}) == 3
                validateattributes(varargin{1}, {'uint16'}, {'size', [nan nan 3]}, mfilename, 'IMDATA');
                data = uint8(varargin{1}*(double(intmax('uint8'))/double(intmax('uint16'))));
            else
                error('imclipboard:InvalideUINT16Image', ...
                    'UINT16 image data must be 2 or 3 dimensions');
            end
            
    end
    
end

end

function filename = validatePasteInput(varargin)

if nargin == 0
    filename = '';
else
    filename = varargin{1};
    validateattributes(filename, {'char'}, {'row'}, mfilename, 'FILENAME');
    [~, ~, e] = fileparts(filename);
    if ~ismember(lower(e), ...
            {'.jpg', '.jpeg', '.tif', '.tiff', '.png', '.gif', '.bmp'})
        error('imclipboard:InvalidFileFormat', ...
            'You must specify a valid image extension: .jpg, .tif, .png, .gif, .bmp');
    end
end

end