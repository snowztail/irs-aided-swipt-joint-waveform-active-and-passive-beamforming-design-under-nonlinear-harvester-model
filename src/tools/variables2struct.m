function s = variables2struct(varargin)
	s = struct;
	for var = 1:nargin
		s.(genvarname(inputname(var), fieldnames(s))) = varargin{var};
	end
end
