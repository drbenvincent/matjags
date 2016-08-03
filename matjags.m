function [samples, stats, structArray] = matjags(dataStruct, jagsModel, initStructs, varargin)
% MATJAGS, a Matlab interface for JAGS
% Version 1.3.1. Tested on JAGS 3.3.0, Windows 64-bit version
%
% This code has been adapted from MATBUGS that was written by Kevin Murphy and Maryam Mahdaviani
%
% [samples, stats] = matjags(dataStruct,  bugsModelFileName,  initStructs, ...)
%
% INPUT:
% dataStruct contains values of observed variables.
% jagsModel is the name of the model file or a string that contains a jags model
% initStructs contains initial values for the latent variables (unlike
% matbugs, this is a required variable)

% Note: variables with names 'a.b' in the bugs model file
% should be called 'a_b' in the matlab data structure.
%
% Optional arguments passed as 'string', value pairs [default in brackets, case
% insensitive]:
% 'monitorParams' - cell array of field names (use 'a_b' instead of 'a.b')
%                   [defaults to *, which currently does nothing...]
% 'nAdapt'   - number of adaptation steps [1000]
% 'nChains'  - number of chains [3]
% 'nBurnin'  - num samples for burn-in per chain [1000]
% 'nSamples' - num samples to keep after burn-in [5000]
% 'thin'     - keep every n'th step [1]
% 'dic'      - [1] read out DIC values
% 'workingDir' - directory to store temporary data/init/coda files [pwd/tmp]. Note that total number of iterations  = nBurnin + nSamples * thin.
% 'savejagsoutput' - 0/1 = do/do not produce text files with output from JAGS
% 'verbosity'  -
%    0 = no messages during runtime;
%    1 = minimum number of messages (e.g. which chain is being executed);
%    2 = maximum number of messages
% 'cleanup' - 0/1 -- do we want to remove PREVIOUS temporary files?
% 'dotranspose' - Set to 0 (default) if you want to insure compatibility with matbugs/WinBUGS
% 'rndseed' - set to 1 to randomise seed for each MCMC chain. Default is
% set to 0, for no randomisation of the seed.
%
% OUTPUT
% S contains the samples; each field may have a different shape:
%  S.theta(c, s)       is the value of theta in sample s, chain c
%                      (scalar variable)
%  S.theta(c, s, i)    is the value of theta(i) in sample s, chain c
%                      (vector variable)
%  S.theta(c, s, i, j) is the value of theta(i,j) in sample s, chain c
%                      (matrix variable)
%
% stats contains various statistics, currently:
%    stats.mean, stats.std and stats.Rhat, stats.DIC.
% Each field may have a different shape:
%    stats.mean.theta
%    stats.mean.theta(i)
%    stats.mean.theta(i,j)
%
% Rhat is the "estimated potential scale reduction" statistic due to
%     Gelman and Rubin.
% Rhat values less than 1.1 mean the chain has probably converged for
%     this variable.
%
% Example
%
% 		[samples, stats ] = matjags( ...
% 			datastruct, ...
% 			fullfile(pwd, 'Gaussian_3.txt'), ...
% 			init0, ...
% 			'doparallel' , 0, ...
% 			'nAdapt', 1000, ...
% 			'nchains', nchains,...
% 			'nburnin', nburnin,...
% 			'nsamples', nsamples, ...
% 			'thin', 1, ...
% 			'dic' , 1,...
% 			'monitorparams', {'mu','sigma'}, ...
% 			'savejagsoutput' , 1 , ...
% 			'verbosity' , 1 , ...
% 			'cleanup' , 0 , ...
% 			'showwarnings' , 1 , ...
% 			'workingdir' , 'tmpjags',...
% 			'rndseed' , 0);
%
% For Windows users:
% The JAGS executable should be placed in the windows path
% In Windows 7, go to Control Panel, System and Security, System
% and click on "Advanced System Settings" followed by "Environment Variables"
% Under System variables, click on Path, and add the jags path to the string
% This could look something like "C:\Program Files\JAGS\JAGS-3.3.0\x64\bin"
%
% For MAC users:
% Mike Kalish provided some suggestions to make matjags work on a Mac. These were implemented
% in the current version but this code is not thouroughly tested yet. Please send me any changes in
% the code needed to make matjags run on a Mac/Linux/Unix system.
%
% Written by Mark Steyvers (mark.steyvers@uci.edu) based on the code
% MATBUGS that was written by Maryam Mahdaviani (maryam@cs.ubc.ca)
% and Kevin Murphy (murphyk@cs.ubc.ca)

% Changes in version 1.3:
% * Bug fix. Added lines marked by "% GP" as suggested by Ganesh Padmanabhan

% Changes in version 1.1:
% * The warnings produces by JAGS are now suppressed by default. This removes
% any message about adaptation being incomplete due to a small number of
% burnin iterations
% * Changed the default working directory for JAGS to make it platform
% independent

defaultworkingDir = tempname();

% Process input options
opts = inputParser;
opts.FunctionName = mfilename;
opts.addRequired('dataStruct',@isstruct);
opts.addRequired('jagsModel',@isstr);
opts.addRequired('initStructs',@isstruct);
opts.addParameter('nChains', 1, @isscalar);
opts.addParameter('workingDir', defaultworkingDir, @isstr);
opts.addParameter('nAdapt',1000, @isscalar);
opts.addParameter('nBurnin',1000, @isscalar)
opts.addParameter('nSamples',5000, @isscalar)
opts.addParameter('monitorParams',{}, @iscellstr)
opts.addParameter('thin',1, @isscalar)
opts.addParameter('dic',1, @isscalar)
opts.addParameter('doParallel',0, @isscalar)
opts.addParameter('savejagsoutput',1, @isscalar)
opts.addParameter('verbosity', 0, @isscalar)
opts.addParameter('cleanup', 0, @isscalar)
opts.addParameter('showwarnings', 0, @isscalar)
opts.addParameter('dotranspose', 0, @isscalar)
opts.addParameter('rndseed', 0, @isscalar)
opts.parse(dataStruct, jagsModel, initStructs, varargin{:});
opts = opts.Results;

% Core matjags logic ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[modelFullPath, workingDirFullPath, isWorkDirTemporary] = set_up();
[jagsDataFullPath, nmonitor] = create_data_file();
make_JAGS_scripts();
[result, status] = run_jags();
error_reporting();
[samples, stats] = coda2matlab();
clean_up();
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	function [modelFullPath, workingDirFullPath, isWorkDirTemporary] = set_up()
		isWorkDirTemporary = strcmp(defaultworkingDir, opts.workingDir) && ~exist(opts.workingDir, 'file');
		
		if length( opts.initStructs ) ~= opts.nChains
			error( 'Number of structures with initial values should match number of chains' );
		end
		
		if is_modelstring(opts.jagsModel)
			workingDirFullPath = get_working_directory(opts.workingDir);
			modelFullPath = fullfile(workingDirFullPath, 'jags_model.jags');
			fid = fopen(modelFullPath, 'w');
			if fid == -1
				error(['Cannot write model to "', modelFullPath, '"' ]);
			end
			fprintf(fid, '%s', opts.jagsModel);
			fclose(fid);
		else
			[modelFullPath, workingDirFullPath] = get_model_and_working_directory_paths(opts.jagsModel, opts.workingDir);
		end
		
		% Do we want to cleanup files before we start?
		if opts.cleanup==1
			delete( fullfile(workingDirFullPath, 'CODA*') );
			delete( fullfile(workingDirFullPath, 'jag*') );
		end
	end

	function [jagsDataFullPath, nmonitor] = create_data_file()
		% Create the data file
		jagsDataFullPath = fullfile(workingDirFullPath, 'jagsdata.R');
		dataGenjags(opts.dataStruct, jagsDataFullPath , '', opts.dotranspose );
		
		nmonitor = length( opts.monitorParams );
		if nmonitor == 0
			error( 'Please specify at least one node name to monitor' );
		end
	end

	function make_JAGS_scripts()
		% Pick a random seed. Remember that 'randi' is itself subject to the random
		% seed, so you may wish to randomise this at the start of your matlab
		% session.
		if opts.rndseed==1
			seed = randi([1 10000000],1);
		end
		
		% Develop a separate JAGS script for each chain
		for whchain=1:opts.nChains
			codastemFullPath     = fullfile(workingDirFullPath, sprintf( 'CODA%d' , whchain ));
			InitDataFullPath     = fullfile(workingDirFullPath, sprintf( 'jagsinit%d.R' , whchain ));
			
			% Create the jags script for this chain
			jagsScriptFullPath   = fullfile(workingDirFullPath, sprintf( 'jagscript%d.cmd' , whchain ));
			[ fid , message ] = fopen( jagsScriptFullPath , 'wt' );
			if fid == -1
				error( message );
			end
			
			if opts.dic
				fprintf( fid , 'load dic\n' );
			end
			
			fprintf( fid , 'model in "%s"\n' , modelFullPath);
			fprintf( fid , 'data in "%s"\n' , jagsDataFullPath );
			fprintf( fid , 'compile, nchains(1)\n' );
			fprintf( fid , 'parameters in "%s"\n' , InitDataFullPath );
			fprintf( fid , 'initialize\n' );
			fprintf( fid , 'adapt %d\n' , opts.nAdapt );
			fprintf( fid , 'update %d\n' , opts.nBurnin );
			for j=1:nmonitor
				fprintf( fid , 'monitor set %s, thin(%d)\n' , opts.monitorParams{ j } , opts.thin );
			end
			if opts.dic
				fprintf( fid , 'monitor deviance\n' );
			end
			fprintf( fid , 'update %d\n' , opts.nSamples * opts.thin );
			fprintf( fid , 'coda *, stem(''%s'')\n' , codastemFullPath );
			fclose( fid );
			
			% Create the init file
			switch opts.rndseed
				case{0}
					addlines = { '".RNG.name" <- "base::Mersenne-Twister"' , ...
						sprintf( '".RNG.seed" <- %d' , whchain ) };
				case{1}
					% Start each chain with a unique random seed.
					addlines = { '".RNG.name" <- "base::Mersenne-Twister"' , ...
						sprintf( '".RNG.seed" <- %d' , whchain+seed ) };
			end
			dataGenjags( opts.initStructs, InitDataFullPath , addlines, opts.dotranspose );
		end
	end

	function [result, status] = run_jags()
			
		status = cell( 1,opts.nChains );
		result = cell( 1,opts.nChains );
			
		% Do we use the Matlab parallel computing toolbox?
		if opts.doParallel
			% open parallel pool
			if isempty(gcp('nocreate'))
				error( 'Matlab pool of workers not initialized. Use command "parpool(7)" for example to open up a pool of 7 workers' );
			end
			parfor whchain=1:opts.nChains
				if opts.verbosity > 0
					fprintf( 'Running chain %d (parallel execution)\n' , whchain  );
				end
				jagsScript   = fullfile(workingDirFullPath, sprintf( 'jagscript%d.cmd' , whchain ));
				[status{ whchain },result{whchain}] = run_jags_script(jagsScript);
			end
		else
			for whchain=1:opts.nChains
				if opts.verbosity > 0
					fprintf( 'Running chain %d (serial execution)\n' , whchain );
				end
				jagsScript   = fullfile(workingDirFullPath, sprintf( 'jagscript%d.cmd' , whchain ));
				[status{ whchain },result{whchain}] = run_jags_script(jagsScript);
			end
		end
		
		% Save the output from JAGS to a text file?
		if opts.savejagsoutput
			for whchain=1:opts.nChains
				filenmFullPath = fullfile(workingDirFullPath, sprintf( 'jagoutput%d.txt' , whchain ));
				[ fid , message ] = fopen( filenmFullPath , 'wt' );
				if fid == -1
					error( message );
				end
				resultnow = result{whchain};
				fprintf( fid , '%s' , resultnow );
				fclose( fid );
			end
		end
	end

	function error_reporting()
		%% Do some error checking.
		% For each chain, check if the output contains some error or warning message.
		for whchain=1:opts.nChains
			resultnow = result{whchain};
			statusnow = status{ whchain };
			if status{whchain} > 0
				error( [ 'Error from system environment: ' resultnow ] );
			end
			
			% Do we get an error message anywhere from JAGS --> produce an error
			pattern = [ 'can''t|RUNTIME ERROR|syntax error|failed' ];
			errstr = regexpi( resultnow , pattern , 'match' );
			if ~isempty( errstr )
				fprintf( 'Error encountered in jags (chain %d). Check output from JAGS below:\n' , whchain  );
				fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
				error( 'Stopping execution because of jags error' );
			end
			
			% Do we get a warning message anywhere from JAGS --> produce a matlab warning
			if opts.showwarnings ~= 0
				pattern = [ 'WARNING' ];
				errstr = regexpi( resultnow , pattern , 'match' );
				if ~isempty( errstr )
					warning( 'JAGS produced a warning message. Check the output below produced by the JAGS run' );
					fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
				end
			end
			
			if opts.verbosity == 2
				fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
			end
			
			% NOTE: if the error is "jags is not recognized as an internal or external
			% command, then the jags bin folder is not on the windows path"
		end
	end

	function [samples, stats] = coda2matlab()
		%% Extract information from the output files so we can pass it back to Matlab
		% the index files are identical across chains, just pick first one
		codaIndexFullPath = fullfile(workingDirFullPath, 'CODA1index.txt');
		for i=1:opts.nChains
			codaFFullPath = fullfile(workingDirFullPath, [ 'CODA' , num2str(i) , 'chain1.txt' ]);
			S = bugs2mat(codaIndexFullPath, codaFFullPath);
			structArray(i) = S;
		end
		samples = structsToArrays(structArray);
		stats = computeStats(samples,opts.dic);
	end

	function clean_up()
		if isWorkDirTemporary
			delete(fullfile(workingDirFullPath, 'jag*'));
			delete(fullfile(workingDirFullPath, 'CODA*'));
			rmdir(workingDirFullPath);
		end
	end

end


%% ----- nested functions -----------------


function result = is_modelstring(string)
result = ~isempty(regexp(string, '^\s*model\s*\{'));
end


function [status, result] = run_jags_script(jagsScript)
if ispc
	jagsPath = 'jags';
else
	possibleDirectories = {'/usr/local/bin/', '/usr/bin/'};
	jagsPath = get_jags_path_from_possible_directories(possibleDirectories);
end
cmd = sprintf('%s %s', jagsPath, jagsScript);
if ispc()
	[status, result] = dos( cmd );
else
	[status, result] = unix( cmd );
end
end


function path = get_jags_path_from_possible_directories(possibleDirectories)
for i=1:length(possibleDirectories)
	if is_jags_directory(possibleDirectories{i})
		path = fullfile(possibleDirectories{i}, 'jags');
		return
	end
end
path = 'jags';
end


function result = is_jags_directory(directory)
if ispc()
	jags = fullfile(directory, 'jags.bat');
else
	jags = fullfile(directory, 'jags');
end
result = exist(jags, 'file');
end


function workingDirFullPath = get_working_directory(workingDir)
curdir = pwd;
% Does the temporary directory exist? If not, create it
if ~exist( workingDir , 'dir' )
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(workingDir);
	if SUCCESS == 0
		error( MESSAGE );
	end
end
cd(workingDir);
workingDirFullPath = pwd();
cd(curdir);
end


function [modelFullPath, workingDirFullPath] = get_model_and_working_directory_paths(jagsFilenm, workingDir)
% get the current directory
curdir = pwd;

[ whdir , jagsModelBase , modelextension ] = fileparts( jagsFilenm );
jagsModel = [ jagsModelBase modelextension ];

cd(whdir);

% expand home dir (~) to absolute path
if strncmp(whdir, '~', 1)
	whdir = [getenv('HOME') whdir(2:end)];
end

if ~isempty(whdir) && (strcmp(whdir(1),filesep) || (length(whdir) > 2 && whdir(2) == ':'))
	% Case when a full path string is specified for the jagsModel
	modelFullPath = fullfile(whdir , jagsModel);
else
	% Case when a relative path string is specified for the jagsModel
	modelFullPath = fullfile(curdir, whdir, jagsModel);
end

% Does the temporary directory exist? If not, create it
if ~exist( workingDir , 'dir' )
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(workingDir);
	if SUCCESS == 0
		error( MESSAGE );
	end
end

cd(workingDir);
workingDirFullPath = pwd();
cd(curdir);
end


function dataGenjags(dataStruct, fileName, addlines, dotranspose )
% This is a helper function to generate data or init files for JAGS
% Inputs:
%   fileName: name of the text file containing initial values. for each
%             chain we'll fileName_i where 'i' is the chain number,
%   dataStruct: is a Struct with name of params(consistant in the same
%               order with paramList) are fields and intial values are functions

if nargin<3
	error(['This function needs three arguments']);
end

fieldNames = fieldnames(dataStruct);
Nparam = size(fieldNames, 1);

fid = fopen(fileName, 'w');
if fid == -1
	error(['Cannot open ', fileName ]);
end

for i=1:Nparam
	fn = fieldNames(i);
	fval = fn{1};
	val = dataStruct.(fval);
	[sfield1, sfield2]= size(val);
	
	msfield = max(sfield1, sfield2);
	newfval = strrep(fval, '_', '.');
	newfval = [ '"' newfval '"' ];
	
	if ((sfield1 == 1) && (sfield2 == 1))  % if the field is a singleton
		fprintf(fid, '%s <-\n%G',newfval, val);
		
		%
		% One-D array:
		%   beta = c(6, 6, ...)
		%
		% 2-D or more:
		%   Y=structure(
		%     .Data = c(1, 2, ...), .Dim = c(30,5))
		%
	elseif ((length(size(val)) == 2) && ((sfield1 == 1) || (sfield2 == 1)))
		fprintf(fid, '%s <-\nc(',newfval);
		for j=1:msfield
			if (isnan(val(j)))
				fprintf(fid,'NA');
			else
				% format for winbugs
				fprintf(fid,wb_strval(val(j)));
			end
			if (j<msfield)
				fprintf(fid, ', ');
			else
				fprintf(fid, ')');
			end
		end
	else
		% non-trivial 2-D or more array
		valsize    = size(val);
		alldatalen = prod(valsize);
		
		%Truccolo-Filho, Wilson <Wilson_Truccolo@brown.edu>
		if length(valsize)<3
			if dotranspose==0
				alldata = reshape(val, [1, alldatalen]);
			else
				alldata = reshape(val', [1, alldatalen]);
				valsize = size( val' );
			end
		elseif length(valsize)==3
			clear valTransp
			if dotranspose==1
				for j=1:valsize(3)
					valTransp(j,:,:)=val(:,:,j)';%need a new variable, since val might be rectangular
				end
				alldata=valTransp(:)';
			else % GP
				alldata = reshape(val, [1, alldatalen]); % GP
			end
		else
			['Error: 4D and higher dimensional arrays not accepted']
			return
		end
		
		fprintf(fid, '%s <-\nstructure(c(', newfval);
		for j=1:alldatalen
			if (isnan(alldata(j)))
				fprintf(fid,'NA');
			else
				% format for winbugs
				fprintf(fid,wb_strval(alldata(j)));
			end
			if (j < alldatalen)
				fprintf(fid,',');
			else
				fprintf(fid,'), .Dim=c(', alldata(j));
			end
		end
		
		for j=1:length(valsize)
			if (j < length(valsize))
				fprintf(fid, '%G,', valsize(j));
			else
				fprintf(fid, '%G))', valsize(j));
			end
		end
	end
	if (i<Nparam)
		fprintf(fid, '\n');
	else
		fprintf(fid, '\n');
	end
end

if length( addlines ) > 0
	nextra = length( addlines );
	for j=1:nextra
		fprintf( fid , '%s\n' , addlines{ j } );
	end
end

fclose(fid);
end


function s = wb_strval(v)
% Converts numeric value to a string that is acceptable by winbugs.
% This is most problematic for exponent values which must have at least 1
% decimal and two digits for the exponent. So 1.0E+01 rather than 1E+001
% Note that only Matlab on PC does 3 digits for exponent.
s = sprintf('%G', v);
if strfind(s, 'E')
	if length(strfind(s, '.')) == 0
		s = strrep(s, 'E', '.0E');
	end
	s = strrep(s, 'E+0', 'E+');
	s = strrep(s, 'E-0', 'E-');
end
end


function A = structsToArrays(S)
% Suppose S is this struct array
%
% S(c).X1(s)
% S(c).X2(s,i)
% S(c).X3(s,i,j)
%
% where s=1:N in all cases
%
% Then we return
% A.X1(c,s)
% A.X2(c,s,i)
% A.X3(c,s,i,j)

C = length(S);
fld = fieldnames(S);
A = [];
for fi=1:length(fld)
	fname = fld{fi};
	tmp = S(1).(fname);
	sz = size(tmp);
	psz = prod(sz);
	data = zeros(C, psz);
	for c=1:C
		tmp = S(c).(fname);
		data(c,:) = tmp(:)';
	end
	if sz(2) > 1 % vector or matrix variable
		data = reshape(data, [C sz]);
	end
	A.(fname) = data;
end
end


function stats = computeStats(all_samples, dodic)
% For each variable (field in the all_samples structure), compute a series
% of statistics (which will be fields of stats). The only complexity is
% that we have to be sensitive to whether each variable is a scalar,
% vector, or 2D matrix. Higher-dimensional variables are not currently
% supported.
disp('Calculating statistics')

variable_names = fieldnames(all_samples);

stats = struct('Rhat',[], 'mean', [], 'median', [], 'std', [],...
	'ci_low' , [] , 'ci_high' , [],...
	'hdi_low', [] , 'hdi_high' , []);

for v=1:length(variable_names)
	var_name = variable_names{v};
	var_samples = all_samples.(var_name);
	
	sz = size(var_samples);
	Nchains = sz(1);
	Nsamples = sz(2);
	dims = ndims(var_samples);
	
	% Calculate stats
	switch dims
		case{2} % scalar
			stats.mode.(var_name)	= calcMode( var_samples(:) );
			stats.median.(var_name) = median( var_samples(:) );
			stats.mean.(var_name)	= mean( var_samples(:) );
			stats.std.(var_name)	= std( var_samples(:) );
			[stats.ci_low.(var_name), stats.ci_high.(var_name)] = calcCI(var_samples(:));
			[stats.hdi_low.(var_name), stats.hdi_high.(var_name)] = calcHDI(var_samples(:));
			
		case{3} % vector
			for n=1:sz(3)
				stats.mode.(var_name)(n)	= calcMode( vec(var_samples(:,:,n)) );
				stats.median.(var_name)(n)	= median( vec(var_samples(:,:,n)) );
				stats.mean.(var_name)(n)	= mean( vec(var_samples(:,:,n)) );
				stats.std.(var_name)(n)		= std( vec(var_samples(:,:,n)) );
				[stats.ci_low.(var_name)(n),...
					stats.ci_high.(var_name)(n)] = calcCI( vec(var_samples(:,:,n)) );
				[stats.hdi_low.(var_name)(n),...
					stats.hdi_high.(var_name)(n)] = calcHDI( vec(var_samples(:,:,n)) );
			end
		case{4} % 2D matrix
			for a=1:sz(3)
				for b=1:sz(4)
					stats.mode.(var_name)(a,b)		= calcMode( vec(var_samples(:,:,a,b)) );
					stats.median.(var_name)(a,b)	= median( vec(var_samples(:,:,a,b)) );
					stats.mean.(var_name)(a,b)		= mean( vec(var_samples(:,:,a,b)) );
					stats.std.(var_name)(a,b)		= std( vec(var_samples(:,:,a,b)) );
					[stats.ci_low.(var_name)(a,b),...
						stats.ci_high.(var_name)(a,b)] = calcCI( vec(var_samples(:,:,a,b)) );
					[stats.hdi_low.(var_name)(a,b),...
						stats.hdi_high.(var_name)(a,b)] = calcHDI( vec(var_samples(:,:,a,b)) );
				end
			end
		otherwise
			warning('calculation of stats not supported for >2D variables. You could implement it and send a pull request.')
			stats.mode.(var_name) = [];
			stats.median.(var_name) = [];
			stats.mean.(var_name) = [];
			stats.std.(var_name) = [];
			stats.ci_low.(var_name) = [];
			stats.ci_high.(var_name) = [];
			stats.hdi_low.(var_name) = [];
			stats.hdi_high.(var_name) = [];
	end
	
	%% "estimated potential scale reduction" statistics due to Gelman and Rubin.
	Rhat = calcRhat();
	if ~isnan(Rhat)
		stats.Rhat.(var_name) = squeeze(Rhat);
	end
	
end

%% DIC calculation
if dodic
	dbar = mean( var_samples.deviance(:));
	dhat = min( var_samples.deviance(:));
	pd = dbar - dhat;
	stats.dic = pd + dbar;
end

	function Rhat = calcRhat()
		st_mean_per_chain = mean(var_samples, 2);
		st_mean_overall   = mean(st_mean_per_chain, 1);
		
		if Nchains > 1
			B = (Nsamples/Nchains-1) * ...
				sum((st_mean_per_chain - repmat(st_mean_overall, [Nchains,1])).^2);
			varPerChain = var(var_samples, 0, 2);
			W = (1/Nchains) * sum(varPerChain);
			vhat = ((Nsamples-1)/Nsamples) * W + (1/Nsamples) * B;
			Rhat = sqrt(vhat./(W+eps));
		else
			warning('Cannot calculate Rhat with 1 chain.')
			Rhat = nan;
		end
	end

	function [low, high] = calcCI(reshaped_samples)
		% get the 95% interval of the posterior
		ci_samples_overall = prctile( reshaped_samples , [ 2.5 97.5 ] , 1 );
		ci_samples_overall_low = ci_samples_overall( 1,: );
		ci_samples_overall_high = ci_samples_overall( 2,: );
		low = squeeze(ci_samples_overall_low);
		high = squeeze(ci_samples_overall_high);
	end

	function [low, high] = 	calcHDI(reshaped_samples)
		% get the 95% highest density intervals of the posterior
		[hdi_samples_overall_low, hdi_samples_overall_high] = HDIofSamples(reshaped_samples);
		low = squeeze(hdi_samples_overall_low);
		high = squeeze(hdi_samples_overall_high);
	end
end


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

data=load(file_out);

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

function mode = calcMode(x)
[F, XI] = ksdensity( x );
[~, ind] = max(F);
mode = XI(ind);
end

function colvec = vec(x)
% turn an input matrix into a column vector
colvec = x(:);
end
