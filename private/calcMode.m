function mode = calcMode(x)
% SLOW MATLAB KSDENSITY
% [F, XI] = ksdensity( x );
% [~, ind] = max(F);
% mode = XI(ind);

% FASTER
if range(x)==0
	% zero variance, so mode = mean. `kde` crashes in this case
	mode = mean(x);
elseif isinf(range(x))
	warning('values supplied seem to have an infitite range')
	mode = [];
else
	[bandwidth,density,xmesh,cdf] = kde(x,...
		1024,...
		min(x)-range(x)/10 ,...
		max(x)+range(x)/10);
	
	[~, index_of_max_value] = max(density);
	mode = xmesh(index_of_max_value);
end
end