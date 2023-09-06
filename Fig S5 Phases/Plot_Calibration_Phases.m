clear all
close all

%% Load phase / coherence CaCO3 + d18O on pmag age model
load('Phase_Coherence_PmagAgeModel.mat');
coherence_080 = 0.336395;
coherence_095 = 0.528167;
coherence_080_plot = [0:0.01*pi():2*pi();ones(1,201).*coherence_080];
coherence_095_plot = [0:0.01*pi():2*pi();ones(1,201).*coherence_095];
obl_phase = -(phase_CaCO3_d18O(phase_CaCO3_d18O == 23.81,2:4) - pi());
obl_coherence = coherence_CaCO3_d18O(phase_CaCO3_d18O == 23.81,2:4);

arch1 = [obl_phase(3):0.001*pi():obl_phase(2);ones(1,110).*obl_coherence(1)];
arch2 = [obl_phase(3):0.001*pi():obl_phase(2);ones(1,110).*obl_coherence(2)];
arch3 = [obl_phase(3):0.001*pi():obl_phase(2);ones(1,110).*obl_coherence(3)];

%% Load and Prepare Data CaCO3 est.
load('CaCO3est.mat');
x = 1:0.1:4;
y_lin = 19.4 .* x -21.114;
y_exp = 5.7525 .* exp(0.6079 .* x);

%% Predefined positions
pos01 = [0.54 0.2 0.4 0.6];
pos02 = [0.07 0.2 0.4 0.6];

%% Generate figure
figure1 = figure;
set(figure1,'colormap',gray);

dcmObj = datacursormode;
set(dcmObj,'UpdateFcn',@custom_cursor);

%Phase CaCO3 ln(Ca/K) - d18O
p1 = polarplot([obl_phase(1) obl_phase(1)],[0 obl_coherence(3)]);
hold all
p2 = polarplot([obl_phase(2) obl_phase(2)],[0 obl_coherence(3)]);
p3 = polarplot([obl_phase(3) obl_phase(3)],[0 obl_coherence(3)]);
p4 = polarplot(arch1(1,:),arch1(2,:));
p5 = polarplot(arch2(1,:),arch2(2,:));
p6 = polarplot(arch3(1,:),arch3(2,:));
p7 = polarplot(coherence_080_plot(1,:),coherence_080_plot(2,:));
p8 = polarplot(coherence_095_plot(1,:),coherence_095_plot(2,:));

set(gca,'ThetaZeroLocation','top',...
    'ThetaTickLabel',{'0',330:-30:30},...
    'position',pos01);
set([p1 p4],'color','k');
set([p5 p6],'color',[0.3 0.3 0.3]);
set([p2 p3],'color',[0.6 0.6 0.6]);
set([p7 p8],'color','r');
set(p7,'linestyle',':');

%ln(Ca/K) - CaCO3
fig(02) = axes('parent',figure1,...
    'position',pos02);
hold(fig(02),'all');
p2(1) = plot(CaCO3_est(1:104,2),CaCO3_est(1:104,1));
p2(2) = plot(CaCO3_est(105:end,2),CaCO3_est(105:end,1));
p2(3) = plot(x,y_lin);
p2(4) = plot(x,y_exp);

%% Figure properties
% set(fig(1),'xlim',[71.2 76.8],...
%     'xgrid','on',...
%     'ylim',[30 70],...
%     'ytick',30:10:60);


%t(1) = title(fig(01),'b) Phase d18O and CaCO3');

set(fig(02),'ylim',[10 60],...
    'xlim',[1 4],...
    'box','on');
set(p2(1:2),'marker','o',...
    'linestyle','none',...
    'markersize',5);
set(p2(1),'color','b',...
    'markeredgecolor','b',...
    'markerfacecolor','b');
set(p2(2),'color','r',...
    'markeredgecolor','r',...
    'markerfacecolor','r');
set(p2(3:4),'color','k');
set(p2(4),'linestyle','--');

yl(2) = xlabel(fig(02),'ln(Ca/K) (van Peer et al., 2017a)');
xl(2) = ylabel(fig(02),'CaCO3 (wt%)');
l(2) = legend(p2(1:4),{'JOIDES','NOCS','Linear fit','Exponential fit'});
t(2) = title(fig(02),'a) Calibration ln(Ca/K) to CaCO3');

set(l(2),'location','northwest');

%% Final figure layout
set(figure1,'units','centimeters',...
    'position',[0 0 21 15],...
    'outerposition',[0 0 21 15],...
    'paperunits','centimeters',...
    'resize','on',...
    'papersize',[21 15],...
    'paperposition',[0 0 21 15],...
    'paperpositionmode','manual');