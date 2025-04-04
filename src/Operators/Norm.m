function sp = Norm(f,type,varargin)
    sp = ScalarProduct(f,f,type,varargin{:});
    sp = sqrt(sp);
end
