function dispim(fn, figh, pos)
% Display bitmap image
%
% dispim(fn, figh, pos)
%
% you will probably get the best result if
% the image is a true colour image (24-bit)
%
% supported formats: jpg, tif, bmp, pcx, hdf, xwd
%
% See also: imread, imwrite, image, axis

x = imread(fn);

if nargin==1
   image(x)
else
   if ~ishandle(figh)
      figure(figh)
   end
   
   if length(pos)==2
      info = imfinfo(fn);
      pos = [pos(1:2) info.Width info.Height];
   end
   ax = axes('parent', figh,...
      'units', 'pixels',...
      'position', pos);
   image(x,...
      'parent', ax)
end

axis image
axis off
