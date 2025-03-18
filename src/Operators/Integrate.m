function int = Integrate(a, varargin)
quadOrder = 2;
if ~isempty(varargin), quadOrder = varargin{1}; end
int = Integrator.compute(a, a.mesh, quadOrder);
end