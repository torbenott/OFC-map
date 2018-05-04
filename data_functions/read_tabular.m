function Data = read_tabular(fname)
%reads .tabluar files from Galazxy into a MATLAB struct with defined
%fieldnames
% fname: file name (full path)
% Data: struct with fieldnames
%       Gene, BaseMean, BaseMeanA, BaseMeanB, foldChange,log2foldChange,
%       pval, padj

if exist(fname,'file')==2
data = readtable(fname,'FileType','text');

%convert to cell array
data = table2cell(data);

%convert fold change and log fold change and p values to numbers (are
%saved in strings) (cols 5-8)

Data.Gene = data(:,1);
Data.BaseMean = cell2mat(data(:,2));
Data.BaseMeanA = cell2mat(data(:,3));
Data.BaseMeanB = cell2mat(data(:,4));
Data.foldChange = str2double(data(:,5));
% Data.log2foldChange = str2double(data(:,6));
Data.pval = str2double(data(:,7));
Data.padj = str2double(data(:,8));

%correct log2foldChange
A = Data.BaseMeanA;
A(A==0)=1;
B = Data.BaseMeanB;
B(B==0)=1;
foldChange = B./A;
log2foldChange = log2(foldChange);
Data.log2foldChange = log2foldChange;

else
    error('read_tabular: file not found.');
end
