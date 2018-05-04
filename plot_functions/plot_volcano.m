function plot_volcano(Data,Params)
%% Plots volcano plot
% Data is struct with expected field names:
%         padj
%         Gene
%         log2foldChange
%         Name  - name of dataset
% params is struct with expected field names:
%         sigp - significant p value cutoff
%         cutp - p value for plot cutoff
%         Axis - axis to plot

yData = -log10(Data.padj(Data.padj<Params.cutp));
xData = Data.log2foldChange(Data.padj<Params.cutp);
significant = Data.padj(Data.padj<Params.cutp) < Params.sigp;
names = Data.Gene(Data.padj<Params.cutp);

scatter(xData,yData,'MarkerFaceColor',[0.2,0.2,0.2],'MarkerEdgeColor','none');
scatter(xData(significant),yData(significant),'MarkerFaceColor',[0.9,0.1,0.1],'MarkerEdgeColor','none');

%axis properties
xlabel('log2(fold change)')
ylabel('-log10(p)')
title(Data.Name)
Params.Axis.XLim = [-(max(abs(xData(abs(xData)~=Inf)))+1), (max(abs(xData(abs(xData)~=Inf)))+1)];
Params.Axis.YLim = [0, (max(abs(yData))+1)];

%label significant genes
sig_idx = find(significant);
for i = 1:length(sig_idx)
    text(Params.Axis.XLim(1),yData(sig_idx(i)),names{sig_idx(i)},'FontSize',6)
end

