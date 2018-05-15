% script to read and parse single cell RNA seq data from GSE92522
% https://www.sciencedirect.com/science/article/pii/S0092867417313119


%where  data file is located
DATAPATH  = 'L:\RetroTrap\Sequencing';
DATAFILE = '1-s2.0-S0092867417313119-mmc1.xlsx';


%read table
XLSData =  readtable(fullfile(DATAPATH,DATAFILE),'ReadRowNames',true);

gene_names = XLSData{2:end,1};
valid = cellfun(@valid_witten_gene,gene_names);

%bring in new format
Data=struct();
Data.Gene = gene_names(valid);
Data.Projection = {'Amyg','NAc','VTA'};
%amygdala
Data.log2foldChange(:,1) = str2double(XLSData{valid,2});
Data.padj(:,1) = str2double(XLSData{valid,3});
%NAc
Data.log2foldChange(:,2) = str2double(XLSData{valid,4});
Data.padj(:,2) = str2double(XLSData{valid,5});
%VTA
Data.log2foldChange(:,3) = str2double(XLSData{valid,6});
Data.padj(:,3) = str2double(XLSData{valid,7});

%our data
Data_IP_file = getDir(fullfile(DATAPATH,'VTA'),'file','all results].tabular');
Data_IP=read_tabular(fullfile(DATAPATH,'VTA',Data_IP_file{1}));

%compare
[intersect_names,idx_witten,idx_ip] = intersect(Data.Gene,Data_IP.Gene,'stable');
ip = Data_IP.log2foldChange(idx_ip);
ip_names = Data_IP.Gene(idx_ip);
ip_log2foldChange = Data_IP.log2foldChange(idx_ip);
wi = Data.log2foldChange(idx_witten,3);
wi_padj = Data.padj(idx_witten,3);
padj_sort = sort(Data_IP.padj(idx_ip));
padj = Data_IP.padj(idx_ip);

%% plot
figure
subplot(1,2,1)
corr_idx = ~isnan(wi) & ~isnan(ip) & abs(ip)<1./eps & abs(wi)<1./eps;
scatter(wi(corr_idx),ip(corr_idx));
[r,p]=corrcoef(wi(corr_idx),ip(corr_idx));
title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])
xlabel('log2 fold change witten data')
ylabel('log2 fold change ip data')

%
ff = find(strcmp('Pou3f1',ip_names));
ip(ff)

% 
ff = find(strcmp('Tnnc1',ip_names));
ip(ff)

%
ff = find(strcmp('Fezf2',ip_names));
ip(ff)

ff = find(strcmp('Mal',ip_names));
ip(ff)


subplot(1,2,2)
corr_idx = ~isnan(wi) & ~isnan(ip) & abs(ip)<1./eps & abs(wi)<1./eps & wi_padj<0.01;
scatter(wi(corr_idx),ip(corr_idx));
[r,p]=corrcoef(wi(corr_idx),ip(corr_idx));
title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])
xlabel('log2 fold change witten data')
ylabel('log2 fold change ip data')

%% functions
function v = valid_witten_gene(name)
v = false;
if ~isempty(name) && ~(strncmp(name,'Gm',2) && ~isnan(str2double(name(3))) ) && isnan(str2double(name(1))) && ~strcmp(name,'genesymbol')
    v=true;
end
end