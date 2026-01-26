function hh = bar3c( varargin )
%BAR3C Extension of bar3, which sets the bar color corresponding to its
%height.
%
% Extra parameter: 'MaxColorVal', maxclr
% This will make the color/height of the bar absolute against this maximum
% value. Otherwise, the colors will be relative against the maximum zdata
% value of all bars on the graph.
%
	
	[abscolor, idxabsc]=getarg('MaxColorVal',varargin{:});
	if idxabsc
		varargin(idxabsc+(0:1))=[];
	end
	
	h = bar3(varargin{:});
	for ii = 1:numel(h)
		zdata = get(h(ii),'Zdata');
		cdata = makecdata(zdata(2:6:end,2),abscolor);
		set(h(ii),'Cdata',cdata, 'facecolor','flat');
	end
	
	if nargout>0, 
		hh = h; 
	end
end

function [val, idx] = getarg(strname,varargin)
	idx = 0;
	val = NaN;
	for jj=1:nargin-2
		if strcmpi(varargin{jj},strname)
			idx = jj;
			val = varargin{jj+1};
			return;
		end
	end
end
function cdata = makecdata(clrs,maxclr)
	cdata = NaN(6*numel(clrs),4);
	for ii=1:numel(clrs)
		cdata((ii-1)*6+(1:6),:)=makesingle_cdata(clrs(ii));
	end
	if nargin>=2
		% it doesn't matter that we put an extra value on cdata(1,1)
		% this vertex is NaN (serves as a separator
		cdata(1,1)=maxclr;
	end
end
function scdata = makesingle_cdata(clr)
	scdata = NaN(6,4);
	scdata(sub2ind(size(scdata),[3,2,2,1,2,4],[2,1,2,2,3,2]))=clr;
end
