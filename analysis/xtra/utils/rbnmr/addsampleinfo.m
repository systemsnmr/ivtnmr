function [A,MSG] = addsampleinfo(A,xlfile,uk)
% ADDSAMPLEINFO: adds info from xl-file to nmr-data struct
%
% SYNTAX:	A = addsampleinfo(A,xlfile,uniquekey)
%
% IN:		A: Cells with data structs as read from rbnmr
%			xlfile: xl-file or {'xlfile', 'sheet'}
%			uniquekey: a key used to identify each sample
%
% OUT:		A: Data struct with the added field 'SampleInfo'
%			MSG: Status message

% The XL-data sheet should have one header row (text) and 
% data describing the samples, with at least one uniqe key:
% The following mappings is valid
% Title => A.Title
% FilePath => A.Info.FilePath
% File => fullfile(A.Info.FilePath, A.Info.FileName)
% DS, DataSet or PlotLabel => A.Info.PlotLabel

% Nils Nyberg, KU, 2008-05-22

%% Define parameter mappings
PM.title	= 'A{i}.Title';
PM.sample	= PM.title;	% Alias for title
PM.label	= PM.title; % Alias for title
PM.filepath = 'A{i}.Info.FilePath';
PM.file		= 'fullfile(A{i}.Info.FilePath,A{i}.Info.FileName)';
PM.path		= 'A{i}.Info.RelativePath';
PM.ds		= 'A{i}.Info.PlotLabel';
PM.dataset	= PM.ds; % Alias for plotlabel
PM.plotlabel = PM.ds;
PM.experiment = PM.ds;

%% Check inputs
if nargin < 2; 
	eval('help addsampleinfo'); error('ADDSAMPLEDATA: Need more inputs'); end
if nargin < 3;
	uk = []; end

if ~iscell(A)
	q = A; clear A; A{1} = q;	
end

%% Check xl-file and read contents
sheet = []; range = [];
if iscell(xlfile)
	try
		tmp = xlfile;
		xlfile = tmp{1};
		sheet = tmp{2};
		range = tmp{3};
	catch ME
		% No action
	end
end

% Open and read xl-file
if isempty([sheet,range])
	try
		[num, txt, raw] = xlsread(xlfile);
	catch ME
		rethrow(ME)
	end
elseif isempty(range)
	try
		[num, txt, raw] = xlsread(xlfile,sheet);
	catch ME
		rethrow(ME)
	end
else
	try
		[num, txt, raw] = xlsread(xlfile,sheet,range);
	catch ME
		rethrow(ME)
	end
end

% find out how many lines to skip in the header
skiplines = 0;
if isnan(raw{skiplines+1,1})
	skiplines = skiplines + 1;
end

% extract the first (non NaN) row as header
XL.header = {raw{skiplines+1,:}};

m = length(XL.header);

% Remove lines with NaN's in 1:st column (empty cells in XL)
n = size(raw,1);
keepthese = ones(n,1);
keepthese(1:skiplines+1) = 0;
for i=skiplines+1+1:n
	if isnan(raw{i,1})
		keepthese(i) = 0;
	end
end

XL.data = raw(keepthese==1,:);

% Transpose the data so it is accessible via header names
for i=1:length(XL.header)
	XLt.(genvarname(XL.header{i})) = XL.data(:,i);
end

%% Find out the unique keys in the file
n = length(XL.data);
uniq_keys = zeros(m,1);
for i=1:m
	try
		uniq_keys(i) = length(unique(XL.data(:,i)))==n;
	catch ME
		% Maybe all numbers?
		try
			uniq_keys(i) = length(unique(cell2mat(XL.data(:,i))))==n;
		catch ME2
			uniq_keys(i) = 0;
		end
	end
end

uniq_keys = find(uniq_keys);
if isempty(uk)
	if isempty(uniq_keys)
		uk = XL.header{1};
		warning(...
			'ADDSAMPLEINFO:noUniqueKey',...
			'No unique headers found, using ''%s'' as key',uk)
	else
		uk = XL.header{uniq_keys(1)};
		disp(sprintf('Testing ''%s'' as key',uk));
	end
end

if ~isfield(XLt,uk)
	error('The ''%s'' header does not seem to exist in the spreadsheet',uk)
end

if ~isempty(uniq_keys) && ~any(strcmp(uk,XL.header(uniq_keys)))
	warning(...
		'ADDSAMPLEINFO:notUniqueField',...
		'The key ''%s'' is not a unique key for all entries in the XL-sheet',...
		uk);
end

%% Determine current parameter mapping
% Check if the unique_key (uk) is a predefined parameter.
% Otherwise, check A.*-/A.Info.*-fields in the first A
if any(strcmpi(uk,fields(PM)));
	PMs = PM.(lower(uk));
elseif any(strcmpi(uk,fields(A{1})))
	PMs = sprintf('A{i}.%s',lower(uk));
elseif any(strcmpi(uk,fields(A{1}.Info)))
	PMs = sprintf('A{i}.Info.%s',lower(uk));
else
	PMs = ' ';
end

% Check the ParameterMapping-string (PMs) is a valid parameter
i=1;
try
	eval(sprintf('ValidParameter = %s;',PMs));
catch ME
	old_uk = uk;
	uk = XL.header{1};
	try
		PMs = PM.(lower(uk));
	catch ME2
		error(...
			'ADDSAMPLEINFO:noSuitableMapping',...
			'No suitable mapping of parameter ''%s'' to any fields in the data.',...
			uk);
	end
	try
		eval(sprintf('ValidParameter = %s;',PMs));
	catch ME2
		warning(...
		'ADDSAMPLEINFO:noParameterMatch',...
		'Parameter ''%s'' is not in the data-cells (i=%d)\nUsing ''%s'' instead',...
		old_uk,i,uk);
	end
	if ~any(strcmp(uk,XL.header(uniq_keys)))
		warning(...
			'ADDSAMPLEINFO:notUniqueField',...
			'The key ''%s'' is not a unique key for all entries in the XL-sheet',...
			uk);
	end
end

if strcmp(' ',PMs) || ~ischar(ValidParameter) || size(ValidParameter,1)>1
	error('This parameter is not valid: ''%s''',PMs);
else
	if nargin < 3
		disp(sprintf('The parameter ''%s'' is mapped to ''%s''',uk,PMs))
	end
end

%% Put the data from the XL-sheet to the correct data cells
% All data is put into the A{i}.SampleInfo-cell

XLtf = fields(XLt);
m = length(XLt.(XLtf{1}));	% No. of datalines in XL-sheet
n = length(A);				% No. of spectra
seen = sparse(n,m);
done = zeros(n,length(XLtf));
str = sprintf('j = strcmp(%s,XLt.%s);',PMs,uk);

% Make sure all entries in the uk-parameter is text, otherwise the 'strcmp'
% command will not work as intended
if ~iscellstr(XLt.(uk))
	for i=1:m
		XLt.(uk){i} = num2str(XLt.(uk){i});
	end
end



for i=1:n
	% Find matching entries	
	j=[];
	try
		eval(str);
	catch ME
		rethrow(ME)
	end
	j = find(j);
	if ~isempty(j);
		
		% Make sure non unique entries are 'distributed' evenly
		seen(i,j) = seen(i,j) + 1;
		ji = find(max(sum(seen(:,j))) == sum(seen(:,j)));
		j = j(ji(1));
		seen(i,j) = seen(i,j) - 1;

		% Add the data, numbers as strings
		for k=1:length(XLtf)
			try
				A{i}.SampleInfo.(XLtf{k}) = num2str(XLt.(XLtf{k}){j});
				done(i,k) = done(i,k) + 1;
			catch ME
				rethrow(ME)
			end
		end
	end
end


	noOfAddedData = sum(done,2);
	noOfAddedData(noOfAddedData==noOfAddedData(1))=1;
	noOfAddedData=sum(noOfAddedData);
	
	noOfAddedParam = sum(done,1);
	noOfAddedParam(noOfAddedParam==noOfAddedParam(1))=1;
	noOfAddedParam=sum(noOfAddedParam);

if nargout == 2;
	MSG.MSG = sprintf(...
		'No of datsets with added info: %d/%d\n', ...
		noOfAddedData,n);
	MSG.MSG = sprintf(...
		'%sNo of added parameters: %d/%d\n',...
		MSG.MSG,noOfAddedParam,length(XLtf));
	MSG.MSG = sprintf(...
		'%sNo of XL-file entries: %dx%d',...
		MSG.MSG,m,length(XLtf));
	MSG.done = done;
end

if noOfAddedData ~= n
	warning(...
		'ADDSAMPLEINFO:notAllDataSets',...
		'Not all datasets were assigned new data from XL-sheet (%d out of %d).',...
		n-noOfAddedData,n);
end


