function apply_group_style(h, g)
	c = {'#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEE'; '#A2142F'};
	l = {'-'; '--'; ':'; '-.'};
	m = {'o'; '+'; 's'; 'x'; '^'; 'v'; '>'; '<'};
	set(h, {'color'}, c(rem(repelem((1 : numel(h) / g) - 1, g), numel(c)) + 1),	{'linestyle'}, l(rem(repmat((1 : g) - 1, [1, numel(h) / g]), numel(l)) + 1), {'marker'}, m(rem((1 : numel(h)) - 1, numel(m)) + 1));
end
