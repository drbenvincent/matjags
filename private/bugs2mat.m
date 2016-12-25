function S=bugs2mat(file_ind,file_out,dir)
%BUGS2MAT  Read (Win)BUGS CODA output to matlab structure
%
% S=bugs2mat(file_ind,file_out,dir)
%  file_ind - index file (in ascii format)
%  file_out - output file (in ascii format)
%  dir      - directory where the files are found (optional)
%  S        - matlab structure, with CODA variables as fields
%
% The samples are stored in added 1'st dimension,
% so that 2 x 3 variable R with 1000 samples would be
% returned as S.R(1000,2,3)
%
% Note1: the data is returned in a structure that makes extraction
% of individual sample sequencies easy: the sequencies are
% directly Nx1 double vectors, as for example S.R(:,1,2).
% The computed statistics must, however, be squeezed,
% as mean(S.R,1) is a 1x2x2 matrix.
%
% Note2: in variable names "." is replaced with "_"

% To change the output structure, edit the 'eval' line in the m-file.
% For example, to return all samples as a cell, wich possibly varying
% number of samples for elements of a multidimensional variable,
% cange the 'eval' line to
%    eval(['S.' varname '={samples};']);
% Then the samples of R(2,1) would be returned as cell S.R(2,1)

% (c) Jouko.Lampinen@hut.fi, 2000
% 2003-01-14 Aki.Vehtari@hut.fi - Replace "." with "_" in variable names
% slightly modified by Maryam Mahdaviani, August 2005 (to suppress redundant output)

if nargin>2
	file_ind=[dir '/' file_ind];
	file_out=[dir '/' file_out];
end

ind=readfile(file_ind);

fprintf('matjags: loading %s... ', file_out)
tic
data=load(file_out);
fprintf('took %3.1f seconds\n', toc)

Nvars=size(ind,1);
S=[];
for k=1:Nvars
	[varname,indexstr]=strtok(ind(k,:));
	varname=strrep(varname,'.','_');
	indices=str2num(indexstr);
	if size(indices)~=[1 2]
		error(['Cannot read line: [' ind(k,:) ']']);
	end
	sdata = size(data);
	%indices
	samples=data(indices(1):indices(2),2);
	varname(varname=='[')='(';
	varname(varname==']')=')';
	leftparen=find(varname=='(');
	outstruct=varname;
	if ~isempty(leftparen)
		outstruct=sprintf('%s(:,%s',varname(1:leftparen-1),varname(leftparen+1:end));
	end
	eval(['S.' outstruct '=samples;']);
end
end


function T=readfile(filename)
f=fopen(filename,'r');
if f==-1, fclose(f); error(filename); end
i=1;
while 1
	clear line;
	line=fgetl(f);
	if ~isstr(line), break, end
	n=length(line);
	T(i,1:n)=line(1:n);
	i=i+1;
end
fclose(f);
end