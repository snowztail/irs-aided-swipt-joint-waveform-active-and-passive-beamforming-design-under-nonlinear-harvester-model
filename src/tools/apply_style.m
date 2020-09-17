function apply_style(h)
	l = {'-'; '--'; ':'; '-.'};
	m = {'o'; '+'; 's'; 'x'; '^'; 'v'};
	set(h, {'linestyle'}, l(rem((1 : numel(h)) - 1, numel(l)) + 1), {'marker'}, m(rem((1 : numel(h)) - 1, numel(m)) + 1));
end
