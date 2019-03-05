function z = dummyFunction(x,varargin)
    z=x;
    for ii = 1:nargin-1
        z = z*varargin{ii};
    end
end