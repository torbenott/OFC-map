%% Make Volcano Plots

%% USER %%%%%%%%%

DATAPATH  = 'L:\RetroTrap\Sequencing';
Subdirs = {'PV','VTA','PAG','VS'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Color',[1,1,1],'Position',[520 ,  252,   869 ,  651]);
nplot = ceil(sqrt(length(Subdirs)));
h_ax = zeros(1,nplot*nplot);

for f = 1:length(Subdirs) 
    
    tabfiles = getDir(fullfile(DATAPATH,Subdirs{f}),'file','all results].tabular');
    
    if ~isempty(tabfiles)
        
        %read data
        Data = read_tabular(fullfile(DATAPATH,Subdirs{f},tabfiles{1}));
        Data.Name = Subdirs{f}; %field expected by plot_volcano()
        
        %create axis
        h_ax = subplot(nplot,nplot,f); hold on
        
        %volcano plot
        Params.cutp = 1;
        Params.sigp = 0.1;
        Params.Axis = h_ax;
            
        plot_volcano(Data,Params)
        
    else
        warning('No data found for %s.',Subdirs{f})
    end
end
RedoTicks(gcf)