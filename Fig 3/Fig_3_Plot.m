clear all
close all

pos01 = [0.10 0.830 0.80 0.16];
pos02 = [0.10 0.680 0.80 0.14];
pos03 = [0.10 0.575 0.80 0.10];
pos04 = [0.10 0.425 0.80 0.14];
pos05 = [0.10 0.320 0.80 0.10];
pos06 = [0.10 0.210 0.80 0.10];
pos07 = [0.10 0.100 0.80 0.10];
pos08 = [0.10 0.070 0.80 0.03];

x_age = [20 28];

interval_1218 = 22.2:0.004:32;
interval_926 = 17.9:0.004:26.5;
interval_1090 = 16:0.004:24.1;
interval_1264 = 17.2:0.004:30.1;
interval_U1406 = 21.8:0.004:26.3;

smoothing = 0.995;

%% La2004
load('La2004_raw.mat');
La2004_age = La2004_ecc(:,1);

% normalised ET curve
La2004_ecc_n = (La2004_ecc(:,2) - mean(La2004_ecc(:,2))) ./ std(La2004_ecc(:,2));
La2004_obl_n = (La2004_obl(:,2) - mean(La2004_obl(:,2))) ./ std(La2004_obl(:,2));
La2004_ET = [La2004_ecc(:,1) La2004_ecc_n + La2004_obl_n];

Data_ecc = timeseries(La2004_ecc(:,2));
Interval_25_05 = [0.002 0.003];
Interval_04_02 = [0.0002 0.0006];
Filter_eccL400_raw = idealfilter(Data_ecc,Interval_25_05,'pass');
Filter_eccL24_raw = idealfilter(Data_ecc,Interval_04_02,'pass');
Filter_eccL400 = [La2004_age(:,1) Filter_eccL400_raw.data ./ std(Filter_eccL400_raw.data)];
Filter_eccL24 = [La2004_age(:,1) Filter_eccL24_raw.data ./ std(Filter_eccL24_raw.data)];

env_obl = [La2004_age mean(La2004_obl(:,2)) + abs(hilbert(La2004_obl(:,2) - mean(La2004_obl(:,2))))];

%% Levy et al 2019 Obliquity sensitivity
load('Levy_etal_2019_Obl_Var');
load('Levy_etal_2019_Ecc_Var');
Levy_etal_2019_Obl_Var(:,1) = Levy_etal_2019_Obl_Var(:,1) ./ 1000;
Levy_etal_2019_Ecc_Var(:,1) = Levy_etal_2019_Ecc_Var(:,1) ./ 1000;

load('Astrochron_d18O_Variances_OM_1000kyrwindow_obl23-27.mat');
Var_U1406_Ecc_1000(:,1) = Var_U1406_Ecc_1000(:,1) ./ 1000;
Var_U1406_Obl_1000(:,1) = Var_U1406_Obl_1000(:,1) ./ 1000;
Var_1090_Ecc_1000(:,1) = Var_1090_Ecc_1000(:,1) ./ 1000;
Var_1090_Obl_1000(:,1) = Var_1090_Obl_1000(:,1) ./ 1000;
Var_1218_Ecc_1000(:,1) = Var_1218_Ecc_1000(:,1) ./ 1000;
Var_1218_Obl_1000(:,1) = Var_1218_Obl_1000(:,1) ./ 1000;
Var_1264_Ecc_1000(:,1) = Var_1264_Ecc_1000(:,1) ./ 1000;
Var_1264_Obl_1000(:,1) = Var_1264_Obl_1000(:,1) ./ 1000;
Var_926_Ecc_1000(:,1) = Var_926_Ecc_1000(:,1) ./ 1000;
Var_926_Obl_1000(:,1) = Var_926_Obl_1000(:,1) ./ 1000;
Var_La2004_Ecc_1000(:,1) = Var_La2004_Ecc_1000(:,1) ./ 1000;
Var_La2004_Obl_1000(:,1) = Var_La2004_Obl_1000(:,1) ./ 1000;

Var_La2004_Ecc_1000(:,2) = Var_La2004_Ecc_1000(:,2) .* 100;
Var_La2004_Obl_1000(:,2) = Var_La2004_Obl_1000(:,2) ./ (40);

%Mask intervals of insuffient sample resolution
Var_1218_Ecc_1000(1:301,1) = NaN;
Var_1218_Obl_1000(1:301,1) = NaN;
Var_1090_Ecc_1000(1001:1551,1) = NaN;
Var_1090_Obl_1000(1001:1551,1) = NaN;
Var_926_Ecc_1000(201:576,1) = NaN;
Var_926_Obl_1000(201:576,1) = NaN;

%% Load isotopes 926
load('Palike_etal_2006_926_Isotopes_CF.mat');
Palike06_926_raw = sortrows([Palike06_926_Age Palike06_926_CF Palike06_926_d13C Palike06_926_d18O Palike06_926_Depth],1);
Palike06_926_raw(:,1) = Palike06_926_raw(:,1)+0.007;
Palike06_926_raw(:,1) = round(Palike06_926_raw(:,1),3);
Palike06_926_agemask = [23.041;25.174;25.665;25.853];
Palike06_926_mask = ismember(Palike06_926_raw(:,1),Palike06_926_agemask);
Palike06_926_raw(Palike06_926_mask) = NaN;

Palike06_926_nan = Palike06_926_raw(~isnan(Palike06_926_raw(:,1)),:);
Palike06_926 = Palike06_926_nan(~isnan(Palike06_926_nan(:,2)),:);
[~,a,~] = unique(Palike06_926(:,1));
Palike06_926 = Palike06_926(a,:);
Palike06_926 = Palike06_926(~isnan(Palike06_926(:,4)),:);

%% Sed Rates
U926_Isotopes_interp_sedrate = [interval_926;interp1(Palike06_926(:,1),Palike06_926(:,5),interval_926,'linear','extrap')]';
[fitobject_U926] = fit(U926_Isotopes_interp_sedrate(:,1),U926_Isotopes_interp_sedrate(:,2),'smoothingspline',...
    'smoothingparam',smoothing);
data_smooth_U926 = [U926_Isotopes_interp_sedrate(:,1) feval(fitobject_U926,U926_Isotopes_interp_sedrate(:,1))];
diff_data_U926 = [(data_smooth_U926(1:end-1,1)+data_smooth_U926(2:end,1))/2 diff(data_smooth_U926(:,2))./0.04];

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

%% Sed Rates
U1218_Isotopes_interp_sedrate = [interval_1218;interp1(Palike06_1218(:,1),Palike06_1218(:,4),interval_1218,'linear','extrap')]';
[fitobject_U1218] = fit(U1218_Isotopes_interp_sedrate(:,1),U1218_Isotopes_interp_sedrate(:,2),'smoothingspline',...
    'smoothingparam',smoothing);
data_smooth_U1218 = [U1218_Isotopes_interp_sedrate(:,1) feval(fitobject_U1218,U1218_Isotopes_interp_sedrate(:,1))];
diff_data_U1218 = [(data_smooth_U1218(1:end-1,1)+data_smooth_U1218(2:end,1))/2 diff(data_smooth_U1218(:,2))./0.04];

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

%% Sed Rates
U1090_Isotopes_interp_sedrate = [interval_1090;interp1(Billups04(:,1),Billups04(:,4),interval_1090,'linear','extrap')]';
[fitobject_U1090] = fit(U1090_Isotopes_interp_sedrate(:,1),U1090_Isotopes_interp_sedrate(:,2),'smoothingspline',...
    'smoothingparam',smoothing);
data_smooth_U1090 = [U1090_Isotopes_interp_sedrate(:,1) feval(fitobject_U1090,U1090_Isotopes_interp_sedrate(:,1))];
diff_data_U1090 = [(data_smooth_U1090(1:end-1,1)+data_smooth_U1090(2:end,1))/2 diff(data_smooth_U1090(:,2))./0.04];

%% 1264
load('Liebrand_etal_2016_1264_Isotopes.mat');
Liebrand16 = sortrows([Liebrand16_Age / 1000 Liebrand16_d13C Liebrand16_d18O Liebrand16_Depth]);
Liebrand16 = Liebrand16(~isnan(Liebrand16(:,2)),:);
[~,a,~] = unique(Liebrand16(:,1));
Liebrand16 = Liebrand16(a,:);
Liebrand16 = Liebrand16(~isnan(Liebrand16(:,3)),:);

%% Sed Rates
U1264_Isotopes_interp_sedrate = [interval_1264;interp1(Liebrand16(:,1),Liebrand16(:,4),interval_1264,'linear','extrap')]';
[fitobject_U1264] = fit(U1264_Isotopes_interp_sedrate(:,1),U1264_Isotopes_interp_sedrate(:,2),'smoothingspline',...
    'smoothingparam',smoothing);
data_smooth_U1264 = [U1264_Isotopes_interp_sedrate(:,1) feval(fitobject_U1264,U1264_Isotopes_interp_sedrate(:,1))];
diff_data_U1264 = [(data_smooth_U1264(1:end-1,1)+data_smooth_U1264(2:end,1))/2 diff(data_smooth_U1264(:,2))./0.04];

%% Load isotopes U1406
Tuning_ETP = xlsread('Tuning_TiePoints.xlsx');
load('U1406_Isotopes.mat');
U1406_Isotopes_Age_raw = [interp1(Tuning_ETP(:,1),Tuning_ETP(:,2),U1406_Isotopes(:,1),'linear','extrap') U1406_Isotopes];

[~,a,~] = unique(U1406_Isotopes_Age_raw(:,1));
U1406_Isotopes_Age = U1406_Isotopes_Age_raw(a,:);

%% Sed Rates
U1406_Isotopes_interp_sedrate = [interval_U1406;interp1(U1406_Isotopes_Age(:,1),U1406_Isotopes_Age(:,2),interval_U1406,'linear','extrap')]';
[fitobject_U1406] = fit(U1406_Isotopes_interp_sedrate(:,1),U1406_Isotopes_interp_sedrate(:,2),'smoothingspline',...
    'smoothingparam',smoothing);
data_smooth_U1406 = [U1406_Isotopes_interp_sedrate(:,1) feval(fitobject_U1406,U1406_Isotopes_interp_sedrate(:,1))];
diff_data_U1406 = [(data_smooth_U1406(1:end-1,1)+data_smooth_U1406(2:end,1))/2 diff(data_smooth_U1406(:,2))./0.04];

%% Sample spacing
range = 16:0.2:30;
for k = 1:length(range)
    index = range(k);
    [id_1090,~] = find(Billups04(:,1)>index);
    if ~isempty(id_1090)
        Billups04_sampling_smooth(k,:) = [index min(id_1090)];
    end
    [id_1218,~] = find(Palike06_1218(:,1)>index);
    if ~isempty(id_1218)
        Palike06_1218_sampling_smooth(k,:) = [index min(id_1218)];
    end
    [id_926,~] = find(Palike06_926(:,1)>index);
    if ~isempty(id_926)
        Palike06_926_sampling_smooth(k,:) = [index min(id_926)];
    end
    [id_1264,~] = find(Liebrand16(:,1)>index);
    if ~isempty(id_1264)
        Liebrand16_sampling_smooth(k,:) = [index min(id_1264)];
    end
    [id_U1406,~] = find(U1406_Isotopes_Age(:,1)>index);
    if ~isempty(id_U1406)
        U1406_sampling_smooth(k,:) = [index min(id_U1406)];
    end
end

sampling_smooth_1090 = [mean([Billups04_sampling_smooth(1:end-1,1),Billups04_sampling_smooth(2:end,1)],2),...
    Billups04_sampling_smooth(2:end,2)-Billups04_sampling_smooth(1:end-1,2)];
sampling_smooth_1218 = [mean([Palike06_1218_sampling_smooth(1:end-1,1),Palike06_1218_sampling_smooth(2:end,1)],2),...
    Palike06_1218_sampling_smooth(2:end,2)-Palike06_1218_sampling_smooth(1:end-1,2)];
sampling_smooth_926 = [mean([Palike06_926_sampling_smooth(1:end-1,1),Palike06_926_sampling_smooth(2:end,1)],2),...
    Palike06_926_sampling_smooth(2:end,2)-Palike06_926_sampling_smooth(1:end-1,2)];
sampling_smooth_1264 = [mean([Liebrand16_sampling_smooth(1:end-1,1),Liebrand16_sampling_smooth(2:end,1)],2),...
    Liebrand16_sampling_smooth(2:end,2)-Liebrand16_sampling_smooth(1:end-1,2)];
sampling_smooth_U1406 = [mean([U1406_sampling_smooth(1:end-1,1),U1406_sampling_smooth(2:end,1)],2),...
    U1406_sampling_smooth(2:end,2)-U1406_sampling_smooth(1:end-1,2)];

sampling_smooth_1218(sampling_smooth_1218(:,1) < 22.1,:) = [];
sampling_smooth_926(sampling_smooth_926(:,1) < 17.8,:) = [];
sampling_smooth_1264(sampling_smooth_1264(:,1) < 17.1,:) = [];
sampling_smooth_U1406(sampling_smooth_U1406(:,1) < 21.8,:) = [];

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
p1(1) = plot(Billups04(:,1),Billups04(:,2),'displayname','1090');
p1(2) = plot(Palike06_1218(:,1),Palike06_1218(:,3),...
    'displayname','1218');
p1(3) = plot(Palike06_926(:,1),Palike06_926(:,4),...
    'displayname','926');
p1(4) = plot(Liebrand16(:,1),Liebrand16(:,3),...
    'displayname','1264');
p1(5) = plot(U1406_Isotopes_Age(:,1),U1406_Isotopes_Age(:,5),...
    'displayname','U1406');

fig(2) = axes('parent',figure1,...
    'position',pos02);
hold(fig(2),'all');
pla2004(1) = plot(Var_La2004_Obl_1000(:,1),Var_La2004_Obl_1000(:,2));
psens(1) = plot(Levy_etal_2019_Obl_Var(:,1),Levy_etal_2019_Obl_Var(:,2));

p2(01) = plot(Var_1090_Obl_1000(:,1),Var_1090_Obl_1000(:,2),...
    'displayname','1090 mid + end');
p2(02) = plot(Var_1218_Obl_1000(:,1),Var_1218_Obl_1000(:,2),...
    'displayname','1218');
p2(03) = plot(Var_926_Obl_1000(:,1),Var_926_Obl_1000(:,2),...
    'displayname','926');
p2(04) = plot(Var_1264_Obl_1000(:,1),Var_1264_Obl_1000(:,2),...
    'displayname','1264');
p2(05) = plot(Var_U1406_Obl_1000(:,1),Var_U1406_Obl_1000(:,2),...
    'displayname','U1406');

fig(4) = axes('parent',figure1,...
    'position',pos03);
hold(fig(4),'all');
p4(1) = plot(La2004_obl(:,1),La2004_obl(:,2));
p4(2) = plot(env_obl(:,1),env_obl(:,2));

fig(3) = axes('parent',figure1,...
    'position',pos04);
hold(fig(3),'all');
pla2004(2) = plot(Var_La2004_Ecc_1000(:,1),Var_La2004_Ecc_1000(:,2));
psens(2) = plot(Levy_etal_2019_Ecc_Var(:,1),Levy_etal_2019_Ecc_Var(:,2));

p3(01) = plot(Var_1090_Ecc_1000(:,1),Var_1090_Ecc_1000(:,2),...
    'displayname','1090 mid + end');
p3(02) = plot(Var_1218_Ecc_1000(:,1),Var_1218_Ecc_1000(:,2),...
    'displayname','1218');
p3(03) = plot(Var_926_Ecc_1000(:,1),Var_926_Ecc_1000(:,2),...
    'displayname','926');
p3(04) = plot(Var_1264_Ecc_1000(:,1),Var_1264_Ecc_1000(:,2),...
    'displayname','1264');
p3(05) = plot(Var_U1406_Ecc_1000(:,1),Var_U1406_Ecc_1000(:,2),...
    'displayname','U1406');

fig(5) = axes('parent',figure1,...
    'position',pos05);
hold(fig(5),'all');
p5(1) = plot(La2004_ecc(:,1),100*La2004_ecc(:,2));
p5(2) = plot(Filter_eccL400(:,1),Filter_eccL400(:,2)+6);
p5(3) = plot(Filter_eccL24(:,1),Filter_eccL24(:,2)./2+9);

figsed(1) = axes('parent',figure1,...
    'position',pos07);
hold(figsed(1),'all');
psite(1) = plot(diff_data_U1090(:,1),diff_data_U1090(:,2));
psite(2) = plot(diff_data_U1218(:,1),diff_data_U1218(:,2));
psite(3) = plot(diff_data_U926(:,1), diff_data_U926(:,2));
psite(4) = plot(diff_data_U1264(:,1),diff_data_U1264(:,2));
psite(5) = plot(diff_data_U1406(:,1),diff_data_U1406(:,2));

figsed(2) = axes('parent',figure1,...
    'position',pos06);
hold(figsed(2),'all');
psite(6) = plot(sampling_smooth_1090(:,1),sampling_smooth_1090(:,2));
psite(7) = plot(sampling_smooth_1218(:,1),sampling_smooth_1218(:,2));
psite(8) = plot(sampling_smooth_926(:,1), sampling_smooth_926(:,2));
psite(9) = plot(sampling_smooth_1264(:,1),sampling_smooth_1264(:,2));
psite(10) = plot(sampling_smooth_U1406(:,1),sampling_smooth_U1406(:,2));
%l(1) = line([18 30],[10 10]); %for every 100 kyr
l(1) = line([18 30],[20 20]); %for every 200 kyr

%plot 8: time scale
fig(6) = axes('parent',figure1,...
    'position',pos08);
hold(fig(6),'all');
pa(01) = patch('Vertices',[min(x_age) min(x_age) 23.03 23.03 ...
    23.03 23.03 max(x_age) max(x_age); 0 1 1 0 0 1 1 0]',...
    'Faces',[1 2 3 4; 5 6 7 8],...
    'FaceColor','none',...
    'edgecolor','k');
te(1) = text(mean([min(x_age),23.03]),0.45,'Miocene');
te(2) = text(mean([max(x_age),23.03]),0.45,'Oligocene');

%% Figure properties
set([fig(1:6) figsed(1:2)],'xlim',[min(x_age) max(x_age)],...
    'xminortick','on',...
    'color','none',...
    'tickdir','out');
fig(6).XAxis.MinorTickValues = min(x_age):0.2:max(x_age);
fig(3).YAxis.MinorTickValues = 0:0.001:0.016;
fig(5).YAxis.MinorTickValues = 0:6;
figsed(2).YAxis.MinorTickValues = 0:20:180;

set([fig(1:5) figsed(1:2)],'xcolor','none',...
    'yminortick','on');
set(fig([4 5]),'yaxislocation','right');

set([fig([2 3]) figsed(2)],'box','on');

set(fig(1),'ylim',[0.4 2.8],...
    'ydir','reverse',...
    'ygrid','on');
set(fig(2),'ylim',[0 0.009],...
    'ytick',0:0.002:0.008,...
    'yticklabel',{0:0.002:0.008},...
    'ygrid','on');
set(fig(3),'ylim',[0 0.018],...
    'ytick',0:0.004:0.016,...
    'ygrid','on');
set(fig(4),'ylim',[22 24.5],...
    'ytick',22:24);
set(fig(5),'ylim',[0 10],...
    'ytick',0:2:6);
set(fig(6),'tickdir','out',...
    'ycolor','none');
% set(figsed,'ylim',[0 110],... %for every 100 kyr
%     'ytick',0:20:100,...
%     'ygrid','on');
set(figsed(2),'ylim',[0 185],... %for every 200 kyr
    'ytick',0:40:160,...
    'ygrid','on',...
    'yaxislocation','right');
set(figsed(1),'ylim',[0 3.5],...
    'ytick',0:3,...
    'ygrid','on');

set(p1(1),'color',[0.6 0.4 1  ]); %1090
set(p1(2),'color',[0   0   1  ]); %1218 
set(p1(3),'color',[0.5 0   0  ]); %926
set(p1(4),'color',[1   0.6 0  ]); %1264
set(p1(5),'color',[1   0.2 0  ]); %U1406

set(p1(1:5),'marker','.',...
    'markersize',4,...
    'linewidth',0.25); 
set(l(1),'color','k');

set([p2(1:5) p3(1:5)],'marker','none');
set([p2(01) p3(01) psite([1 6])],'color',get(p1(01),'color'),...
    'linestyle','-');
set([p2(02) p3(02) psite([2 7])],'color',get(p1(02),'color'),...
    'linestyle','-');
set([p2(03) p3(03) psite([3 8])],'color',get(p1(03),'color'),...
    'linestyle','-');
set([p2(04) p3(04) psite([4 9])],'color',get(p1(04),'color'),...
    'linestyle','-');
set([p2(05) p3(05) psite([5 10])],'color',get(p1(05),'color'),...
    'linestyle','-');

set(pla2004(1:2),'color','k',...
    'linestyle','-.');

set(psite(6:10),'marker','.',...
    'markersize',10);

%Levy et al variances
set(psens(1:2),'color','k','linestyle','-');

set(p4(1),'color',[0.8500, 0.3250, 0.0980]);
set(p5(1),'color',[0.3010, 0.7450, 0.9330]);
set([p4(2) p5(2)],'color',[0.3 0.3 0.3]);
set(p5(3),'color','k','linestyle','--');

yl(1) = ylabel(figsed(1),{'Sed. Rate','(cm/kyr)'});
yl(2) = ylabel(figsed(2),{'Sample Resolution','(samples/200 kyr)'});
yl(3) = ylabel(fig(1),'d18O');
yl(4) = ylabel(fig(2),{'Obl. Var. (permille2)'});
yl(5) = ylabel(fig(3),{'Ecc. Var. (permille2)'});
yl(6) = ylabel(fig(4),{'Obliquity','(degrees)'});
yl(7) = ylabel(fig(5),{'Eccentricity','(%)'});


xl(1) = xlabel(fig(6),'Age (Myr ago)');

legend([p2([5 3 4 1 2]) psens(1) pla2004(1)],...
    {'U1406 N Atl','926 Eq Atl','1264 S Atl','1090 S Atl','1218 Eq Pac'...
    'L19 / DV17','La2004,Vobl'});

%% Final layout
set(figure1,'units','centimeters',...
    'position',[0 0 20 24],...
    'outerposition',[0 0 20 24],...
    'paperunits','centimeters',...
    'resize','on');