%% Author(s)

% Corentin Wicht (script, protocol) 
% GitHub : https://github.com/CorentinWicht

% If you have questions or want to contribute to this pipeline, feel free 
% to contact corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% Function to plot t-test results

function Plot_test(EEG,CompareFN,CurrentDesign,Stats,Pval,TempData,Normalization)

% Thresholding original data with P-val (user defined)
TempP = Stats.Original.P;
TempP(TempP > Pval) = 0;

% FDR Correction (should it really only be applied to
% significant ones??)
[~,pN] = limo_FDR(TempP,Pval);

% PLOT TEST %%%%%%%%%%%%%%
h =  findobj('type','figure');
figure(length(h)+1); clf
annotation('textbox',[.35 .5 .5 .5],'String',...
strrep(Stats.FileName,'_','-'),'FitBoxToText','on','FontSize',20,'EdgeColor','w');
Labels = {EEG.chanlocs.labels};
Dat = Stats.Original.P;
LogicalDat = TempP;
LogicalDat(TempP>0)=1;
XTStr = round(EEG.xmin*1000):50:round(EEG.xmax*1000,-1);
Ticks = linspace(1,size(Dat,2),size(XTStr,2));
MinVal=min([min(min(min(TempData.(CompareFN{1}))))...
    min(min(min(TempData.(CompareFN{2}))))]);
MaxVal=max([max(max(max(TempData.(CompareFN{1}))))...
    max(max(max(TempData.(CompareFN{2}))))]);
for pp=1:2
    subplot(2,2,pp)
    PlotData = squeeze(mean(permute(TempData.(CompareFN{pp}),[2 3 1]),3));
    imagesc(PlotData)
    title(CurrentDesign{pp+2})
    set(gca,'ydir','norm')
    hcb = colorbar;
    if strcmpi(Normalization,'Y')
        set(get(hcb,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
    else
        set(get(hcb,'XLabel'),'String','Power (\muV)');
    end
    caxis([MinVal MaxVal])
    colormap(gca,'hsv(50)')
    set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
    set(gca,'YTick',pp:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(pp:2:end));
    hold on; contour(logical(LogicalDat),1,'linecolor','k'); hold off;
end

subplot(2,2,3)
h = imagesc(TempP);
set(gca,'ydir','norm')
title('Stats (p-val)')
colorbar
colormap(gca,'hot(50)')
set(h,'alphadata',TempP~=0)
set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

subplot(2,2,4)
DatFDR=Dat;
DatFDR(DatFDR>pN)=0;
hh=imagesc(DatFDR);
set(gca,'ydir','norm')
title('FDR-corrected P-val')
colorbar
colormap(gca,'hot(50)')
set(hh,'alphadata',DatFDR~=0)
set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

end