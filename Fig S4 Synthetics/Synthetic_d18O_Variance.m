clear all
close all

pos01 = [0.10 0.830 0.80 0.16];
pos02 = [0.10 0.680 0.80 0.14];
pos03 = [0.10 0.575 0.80 0.10];
pos04 = [0.10 0.425 0.80 0.14];

x_age = [22 26];
interval_1218 = 22.2:0.004:32;
interval_1090 = 16:0.004:24.1;

smoothing = 0.995;

%% Levy et al 2019 Obliquity sensitivity
load('Levy_etal_2019_Obl_Var');
Levy_etal_2019_Obl_Var(:,1) = Levy_etal_2019_Obl_Var(:,1) ./ 1000;

load('Astrochron_d18O_Variances_OM_1000kyrwindow_obl23-27.mat');
Var_1090_Ecc_1000(:,1) = Var_1090_Ecc_1000(:,1) ./ 1000;
Var_1090_Obl_1000(:,1) = Var_1090_Obl_1000(:,1) ./ 1000;
Var_1218_Ecc_1000(:,1) = Var_1218_Ecc_1000(:,1) ./ 1000;
Var_1218_Obl_1000(:,1) = Var_1218_Obl_1000(:,1) ./ 1000;

%Mask intervals of insuffient sample resolution
Var_1218_Obl_1000(1:301,1) = NaN;
Var_1090_Obl_1000(1001:1551,1) = NaN;

%% 1218
load('Palike_etal_2006_1218_Isotopes.mat');
Palike06_1218_Depth_raw = [Palike06_1218_Depth; Palike06_1218_Depth; Palike06_1218_Depth];
Palike06_1218_Age_raw = [Palike06_1218_Age_P06m/1000; Palike06_1218_Age_P06m/1000; Palike06_1218_Age_P06m/1000];
Palike06_1218_d18O_raw = [Palike06_1218_d18O_Cspp; Palike06_1218_d18O_Chav; Palike06_1218_d18O_Cgri];
Palike06_1218_d13C_raw = [Palike06_1218_d13C_Cspp; Palike06_1218_d13C_Chav; Palike06_1218_d13C_Cgri];
Palike06_1218_raw = sortrows([Palike06_1218_Age_raw Palike06_1218_d13C_raw Palike06_1218_d18O_raw Palike06_1218_Depth_raw],1);

Palike06_1218_raw(Palike06_1218_raw(:,1) < 22.18,1) = NaN;
Palike06_1218_agemask = [22.241;22.801;23.034;23.523;25.284;25.928;27.735];
Palike06_1218_mask = ismember(Palike06_1218_raw(:,1),Palike06_1218_agemask);
Palike06_1218_raw(Palike06_1218_mask) = NaN;

Palike06_1218_nan = Palike06_1218_raw(~isnan(Palike06_1218_raw(:,1)),:);
Palike06_1218 = Palike06_1218_nan(~isnan(Palike06_1218_nan(:,2)),:);

[~,a,~] = unique(Palike06_1218(:,1));
Palike06_1218 = Palike06_1218(a,:);
Palike06_1218 = Palike06_1218(~isnan(Palike06_1218(:,2)),:);

%% 1090
load('Billups_etal_2004_1090_Isotopes.mat');
Billups04_Age = round(Billups04_Age + 0.0072,4);
Billups04_raw = sortrows([Billups04_Age Billups04_d18O Billups04_d13C Billups04_Depth],1);
Billups04_agemask = [20.0852;20.6092;20.6192;20.6292;20.7272;20.8122;20.8582;20.9122;20.9562;21.0392;21.1132;21.1762;21.2102;21.2152;21.7512;24.1392];
Billups04_mask = ismember(Billups04_raw(:,1),Billups04_agemask);
Billups04_raw(Billups04_mask) = NaN;
Billups04_nan = Billups04_raw(~isnan(Billups04_raw(:,1)),:);
Billups04 = Billups04_nan(~isnan(Billups04_nan(:,2)),:);

[~,a,~] = unique(Billups04(:,1));
Billups04 = Billups04(a,:);
Billups04 = Billups04(~isnan(Billups04(:,2)),:);

%% Sine waves
sine05 = [22:0.001:26;0.5*sin((22:0.001:26).*2*pi()*(1/0.041))]';
sine1 = [22:0.001:26;sin((22:0.001:26).*2*pi()*(1/0.041))]';
sine_combined = [sine1(1:2000,:);sine05(2001:end,:)];

%% Synthetic variances
[Synthetics_upload,~,~] = xlsread('SI_Synthetic_d18O_Variance.xlsx');
synthetics1 = Synthetics_upload(:,1:2);
synthetics05 = Synthetics_upload(:,6:7);
synthetics_comb = Synthetics_upload(:,11:12);

%% Figure
figure1 = figure('RendererMode','manual',...
    'Renderer','Painters',...
    'InvertHardCopy','off');
dcmObj = datacursormode;
set(dcmObj,'UpdateFcn',@custom_cursor);
colormap('jet');

fig(1) = axes('parent',figure1,...
    'position',pos01);
hold(fig(1),'all');
p1(1) = plot(sine1(:,1),sine1(:,2));
p1(2) = plot(sine05(:,1),sine05(:,2));
%p1(3) = plot(sine_combined(:,1),sine_combined(:,2));

fig(2) = axes('parent',figure1,...
    'position',pos02);
hold(fig(2),'all');
p2(1) = plot(synthetics1(:,1)/1000,synthetics1(:,2));
p2(2) = plot(synthetics05(:,1)/1000,synthetics05(:,2));
p2(3) = plot(synthetics_comb(:,1)/1000,synthetics_comb(:,2));

fig(3) = axes('parent',figure1,...
    'position',pos03);
hold(fig(3),'all');
p3(03) = plot(Levy_etal_2019_Obl_Var(:,1),Levy_etal_2019_Obl_Var(:,2));
p3(01) = plot(Var_1090_Obl_1000(:,1),Var_1090_Obl_1000(:,2),...
    'displayname','1090');
p3(02) = plot(Var_1218_Obl_1000(:,1),Var_1218_Obl_1000(:,2),...
    'displayname','1218');

fig(4) = axes('parent',figure1,...
    'position',pos04);
hold(fig(4),'all');
p4(1) = plot(Billups04(:,1),Billups04(:,2),'displayname','1090');
p4(2) = plot(Palike06_1218(:,1),Palike06_1218(:,3),...
    'displayname','1218');


%% Figure properties
set(fig(1:4),'xlim',[min(x_age) max(x_age)],...
    'xminortick','on',...
    'color','none',...
    'tickdir','out',...
    'yminortick','on');

set(fig(1:3),'xcolor','none');
set(fig([2 4]),'yaxislocation','right');

set(fig(1),'ylim',[-1.1 1.1]);
set(fig(2),'ylim',[0 0.6]);
set(fig(3),'ylim',[0 0.01]);
set(fig(4),'ylim',[0.7 2.7],...
    'ydir','reverse',...
    'ytick',1:0.5:2.5);

set(p3(1),'color',[0.6 0.4  1  ]); %1090
set(p3(2),'color',[0.1 0.86 1  ]); %1218 
set(p3(3),'color','k','linestyle','-');

set(p1(1),'color',get(p3(1),'color'));
set(p1(2),'color',get(p3(2),'color'));
set(p2(1),'color',get(p3(1),'color'));
set(p2(2),'color',get(p3(2),'color'));
set(p2(3),'color',get(p3(3),'color'));
set(p4(1),'color',get(p3(1),'color'));
set(p4(2),'color',get(p3(2),'color'));

yl(1) = ylabel(fig(1),{'Synthetic','Sine Waves (-)'});
yl(2) = ylabel(fig(2),{'Obl. Var. Synthetics (-)'});
yl(3) = ylabel(fig(3),{'Obl. Var. 1090/1218 (permille2)'});
yl(4) = ylabel(fig(4),{'d18Ob permille'});

xl(1) = xlabel(fig(4),'Age (Myr ago)');

l(1) = legend(fig(1),{'0.5*Sin','1.0*Sin'});
l(2) = legend(fig(2),{'0.5*Sin','1.0*Sin','Combined'});
l(3) = legend(fig(3),{'Levy et al. (2019)','1090','1218'});
l(4) = legend(fig(4),{'1090','1218'});

%% Final layout
set(figure1,'units','centimeters',...
    'position',[0 0 20 24],...
    'outerposition',[0 0 20 24],...
    'paperunits','centimeters',...
    'resize','on');