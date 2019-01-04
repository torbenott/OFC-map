
%%%%%%DATA
BASE = 'C:\Data\DataPostdoc\MAPseq\data_ZL161';
datafile = 'data_ZL161.mat'; %Longwen's mat file
xlsfile = 'SampleList.xlsx'; %included in repo

%%%%%%%PARAMS
thr=2; %Molecule count threshold (has to be greater than this value)
thr_barcode = 2; %number of unique barcodes per area threshold (has to be greater than this value)

%MAPseq data
data=load(fullfile(BASE,datafile));

%xls sheet for sample info
[~,c]=xlsread(xlsfile,'sample information','H2:H119');
d=c(cellfun(@(x) ~isempty(x),c));

%areas
areas = unique(d);
area_idx=cellfun(@(x) find(ismember(areas,x)),d);
ofc_idx = find(strcmp(areas,'OFC'));

%grand sum (barcodes in OFC AND any other target area
valid = any(data.barcodematrix(:,area_idx==ofc_idx)>thr,2) & any(data.barcodematrix(:,area_idx~=ofc_idx)>thr,2);
NBARCODE = sum(valid);

%orphans
for t = 0:10
orphan = ~any(data.barcodematrix(:,area_idx==ofc_idx)>t,2) & any(data.barcodematrix(:,area_idx~=ofc_idx)>t,2);
NORPHAN=sum(orphan);
PORPHAN(t+1) = NORPHAN./size(data.barcodematrix,1);
end
figure,plot(0:10,PORPHAN,'-.k')
xlabel('Threshold (molecule count)');ylabel('% orphan')

%negative control
neg_idx = find(strcmp(areas,'N'));
n_neg = sum(any(data.barcodematrix(:,area_idx==neg_idx)>0,2));
count_neg = data.barcodematrix(any(data.barcodematrix(:,area_idx==neg_idx)>0,2),area_idx==neg_idx);
fprintf('Negative control shows %i unique barcodes with following counts:\n',n_neg);
disp(count_neg)


%% pairwise overlap over/under representation P(A,B)/P(A)
PP=[];
for i = 1 : length(areas)
    iidx=any(data.barcodematrix(:,area_idx==i),2)&valid;
    for k = 1 : length(areas)
        kidx=any(data.barcodematrix(:,area_idx==k)>thr,2)&valid;
        PP(i,k)=sum(iidx&kidx)./NBARCODE; %P(A,B)
    end
end

MM = PP./repmat(diag(PP),1,size(PP,1)); %P(A,B)/P(A)

%exclude source area and low unique barcode numbers
include = NAREA>thr_barcode & NAREA<max(NAREA);
areas_clean = areas(include);%exlucde zeros and source area
MM_clean = MM(include,include);

%plot
figure
ax=imagesc(MM_clean);
ax.Parent.YTick=1:length(areas_clean);
ax.Parent.YTickLabel = areas_clean;
ax.Parent.XTick=1:length(areas_clean);
ax.Parent.XTickLabel = areas_clean;
colorbar
title('P(A,B)/P(A)')
xlabel('B');ylabel('A')
