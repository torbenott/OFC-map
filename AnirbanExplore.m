% script to read and parse single cell RNA seq data from GSE92522
% https://www.cell.com/cell/comments/S0092-8674(17)30990-X
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92522

%where Anirban's data file is located
DATAPATH  = 'L:\RetroTrap\Sequencing';
DATAFILE = 'GSE92522_paul2016.processed.data.txt';

%read table
Data =  readtable(fullfile(DATAPATH,DATAFILE),'FileType','text','ReadRowNames',true);

sample_names = Data.Properties.VariableNames;
gene_names = Data.Properties.RowNames;

%get pv samples
pv_idx = cellfun(@(x) strncmp(x,'PVC',3), sample_names);

%calculate log2 fold change and p-value (not adjusted) for PV samples for
%each gene
p_pv=nan(1,length(gene_names));
fChange = nan(1,length(gene_names));
parfor i =1 : length(gene_names)
    d = Data{i,:};
    p_pv(i) = ranksum(d(pv_idx),d(~pv_idx));
    fChange(i) = log2(mean(d(pv_idx))./mean(d(~pv_idx)));
end



% load our pv seq data
pv_files = getDir(fullfile(DATAPATH,'PV'),'file','all results].tabular');
pv_data=read_tabular(fullfile(DATAPATH,'PV',pv_files{1}));

%perspective of anirban's data
%sort by significance
[~, sortidx]=sort(p_pv);
gene_names_sort = gene_names(sortidx);
% gene_names_sort(1:20)
n_best = 20;
names_comp = gene_names_sort(1:n_best);
gene_idx = cellfun(@(x) find(strcmp(x,pv_data.Gene)),names_comp,'UniformOutput',false);
% gene_idx(cellfun(@isempty,gene_idx)) = {NaN};
gene_idx=cell2mat(gene_idx);

for i =1:length(gene_idx)
    fprintf('%s %f %f\n',pv_data.Gene{gene_idx(i)}, pv_data.log2foldChange(gene_idx(i)), pv_data.padj(gene_idx(i)))
end

%perspective of our data
[p_pv_sort, sortidx]=sort(pv_data.padj);
gene_names_sort = pv_data.Gene(sortidx);
% gene_names_sort(1:20)
n_best = 30;
names_comp = gene_names_sort(1:n_best);
gene_idx = cellfun(@(x) find(strcmp(x,gene_names)),names_comp,'UniformOutput',false);
% gene_idx(cellfun(@isempty,gene_idx)) = {NaN};
gene_idx=cell2mat(gene_idx);

for i =1:length(gene_idx)
    fprintf('%s %f %f\n',gene_names{gene_idx(i)}, fChange(gene_idx(i)), p_pv(gene_idx(i)))
end

%correlation of all genes
gene_idx_all = cellfun(@(x) find(strcmp(x,pv_data.Gene)),gene_names,'UniformOutput',false);
gene_idx_all(cellfun(@isempty,gene_idx_all)) = {NaN};
gene_idx_all=cell2mat(gene_idx_all);

sc = fChange(~isnan(gene_idx_all)); sc = sc(:);
ip = pv_data.log2foldChange(gene_idx_all(~isnan(gene_idx_all))); ip = ip(:);
ip_p = pv_data.padj(gene_idx_all(~isnan(gene_idx_all))); ip_p = ip_p(:);

figure
subplot(1,2,1)
%exclude genes?
nan_idx = isnan(sc) | isnan(ip) | abs(sc)==Inf | abs(ip) == Inf | abs(ip)<eps | abs(sc)<eps | ip_p>p_pv_sort(50);
scatter(sc(~nan_idx),ip(~nan_idx))
[r,p] = corrcoef(sc(~nan_idx),ip(~nan_idx));
xlabel('log2 fold change single cell data')
ylabel('log2 fold change ip data')
title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])

subplot(1,2,2)
%exclude genes?
nan_idx = isnan(sc) | isnan(ip) | abs(sc)==Inf | abs(ip) == Inf | abs(ip)<eps | abs(sc)<eps;
scatter(sc(~nan_idx),ip(~nan_idx))
[r,p] = corrcoef(sc(~nan_idx),ip(~nan_idx));
xlabel('log2 fold change single cell data')
ylabel('log2 fold change ip data')
title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])
