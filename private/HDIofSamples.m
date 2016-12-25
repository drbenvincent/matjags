function [HDI_lower, HDI_upper] = HDIofSamples(samples)
% Calculate the 95% Highest Density Intervals. This has advantages over the
% regular 95% credible interval for some 'shapes' of distribution.
%
% Translated by Benjamin T. Vincent (www.inferenceLab.com) from code in:
% Kruschke, J. K. (2015). Doing Bayesian Data Analysis: A Tutorial with R,
% JAGS, and Stan. Academic Press.

credibilityMass = 0.95;

[nSamples, N] = size(samples);
for i=1:N
	selectedSortedSamples = sort(samples(:,i));
	ciIdxInc = floor( credibilityMass * numel( selectedSortedSamples ) );
	nCIs = numel( selectedSortedSamples ) - ciIdxInc;
	
	ciWidth=zeros(nCIs,1);
	for n =1:nCIs
		ciWidth(n) = selectedSortedSamples( n + ciIdxInc ) - selectedSortedSamples(n);
	end
	
	[~, minInd] = min(ciWidth);
	HDI_lower(i)	= selectedSortedSamples( minInd );
	HDI_upper(i)	= selectedSortedSamples( minInd + ciIdxInc);
end
end