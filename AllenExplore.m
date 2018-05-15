%allen frontal genes layer V

%where  data file is located
DATAPATH  = 'L:\RetroTrap\Sequencing';
DATAFILENEG = 'AllenFrontalGenesNeg.xml'; %depleted in layer V
DATAFILEPOS = 'AllenFrontalGenesPos.xml'; %enriched in layer V

XMLDataNeg = xml2struct(fullfile(DATAPATH,DATAFILENEG));
XMLDataPos = xml2struct(fullfile(DATAPATH,DATAFILEPOS));

%depleted layer V genes
n1 = numel( XMLDataNeg.Response.objects.object);

Data=struct('Gene',[],'log2foldChange',[]);
for i = 1:n1
    name = XMLDataNeg.Response.objects.object{i}.gene_dash_symbol.Text;
    fold =str2double(XMLDataNeg.Response.objects.object{i}.fold_dash_change.Text);
    Data.Gene = [Data.Gene; {name}];
    Data.log2foldChange=[Data.log2foldChange;-fold];
end

%enriched layer V genes
n2 = numel( XMLDataPos.Response.objects.object);

% Data=struct('Gene',[],'log2foldChange',nan(n2,1));
for i = 1:n2
    name = XMLDataPos.Response.objects.object{i}.gene_dash_symbol.Text;
    fold =str2double(XMLDataPos.Response.objects.object{i}.fold_dash_change.Text);
    Data.Gene = [Data.Gene; {name}];
    Data.log2foldChange=[Data.log2foldChange;fold];
end

%% our data
Subdirs = {'PAG'};
Data_IP=[];
for i =1:length(Subdirs)
Data_IP_file = getDir(fullfile(DATAPATH,Subdirs{i}),'file','all results].tabular');
Data_IP_single=read_tabular(fullfile(DATAPATH,Subdirs{i},Data_IP_file{1}));
Data_IP = [Data_IP,Data_IP_single];
end

%compare
[intersect_names,idx_allen,idx_ip] = intersect(Data.Gene,Data_IP(1).Gene,'stable');
ip = Data_IP(1).log2foldChange(idx_ip);
% ip2 = Data_IP(2).log2foldChange(idx_ip);
% ip3 = Data_IP(3).log2foldChange(idx_ip);
% ip = mean([ip1,ip2,ip3],2);
ip_names = Data_IP(1).Gene(idx_ip);

all = Data.log2foldChange(idx_allen);

% padj_sort = sort(Data_IP.padj(idx_ip));
% padj = Data_IP.padj(idx_ip);

%% plot
figure
subplot(1,2,1)
corr_idx = ~isnan(all) & ~isnan(ip);
scatter(all(corr_idx),ip(corr_idx));
[r,p]=corrcoef(all(corr_idx),ip(corr_idx));
title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])
xlabel('fold change allen data')
ylabel('log2 fold change ip data')
axis square

name = ip_names(ip>5);
ff = find(strcmp(name{1},Data_IP(1).Gene));
Data_IP(1).log2foldChange(ff)
Data_IP(1).padj(ff)


% subplot(1,2,2)
% corr_idx = ~isnan(all) & ~isnan(ip) & abs(ip)<1./eps & abs(all)<1./eps & all_padj<0.01;
% scatter(all(corr_idx),ip(corr_idx));
% [r,p]=corrcoef(all(corr_idx),ip(corr_idx));
% title(['r = ',num2str(r(1,2)),', p = ',num2str(p(1,2))])
% xlabel('log2 fold change witten data')
% ylabel('log2 fold change ip data')