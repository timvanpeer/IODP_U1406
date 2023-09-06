%clear all
%close all

original_path = cd();

pos01 = [0.120 0.080 0.805 0.140]; % U1406 ln(Ca/K)
pos02 = [0.120 0.150 0.805 0.140]; % U1406 d18O
pos03 = [0.120 0.250 0.805 0.140]; % Filtered 18-28 Myr-1 on ln(Ca/K) & obl
pos04 = [0.120 0.340 0.805 0.140]; % Amplitudes obliquity filters
pos05 = [0.120 0.420 0.805 0.140]; % Filtered 6-12 Myr-1 on d18O & ecc
pos06 = [0.120 0.515 0.805 0.140]; % Filtered 2-3 Myr-1 on d18O & ecc
pos07 = [0.120 0.660 0.805 0.165]; % Wavelet ln(Ca/K) U1406
pos08 = [0.925 pos07(2) 0.070 pos07(4)]; % Wavelet Global Power Spectrum ln(Ca/K)
pos09 = [0.120 0.830 0.805 0.165]; % Wavelet d18O U1406
pos10 = [0.925 pos09(2) 0.070 pos09(4)]; % Wavelet Global Power Spectrum d18O

xpos_age = [21.8 22:0.5:26 26.2];

%% Load tuning spreadsheet
Tuning_ETP = xlsread('Tuning_TiePoints.xlsx');

%% Load and Prepare ln(Ca/K) Data U1406
lnCaK_depth = xlsread('Data_lnCaK_depth.xlsx');

CaCO3_depth = [lnCaK_depth(:,1) 5.7525 .* exp(0.6079 .* lnCaK_depth(:,2))];

lnCaK_ETP_raw = [interp1(Tuning_ETP(:,1),Tuning_ETP(:,2),CaCO3_depth(:,1),'linear','extrap') CaCO3_depth(:,2)];
lnCaK_ETP = [21.04:0.001:27.1;interp1(lnCaK_ETP_raw(:,1),lnCaK_ETP_raw(:,2),21.04:0.001:27.1,'linear','extrap')]';
Data_notchfilter_lnCaK_ETP = timeseries(lnCaK_ETP(:,2));
Interval_notch_ETP = [0 0.001];
Notchfilter_lnCaK_ETP = idealfilter(Data_notchfilter_lnCaK_ETP,Interval_notch_ETP,'notch');
Notch_lnCaK_ETP = [lnCaK_ETP(:,1) Notchfilter_lnCaK_ETP.data];
Notch_lnCaK_ETP_plot = [lnCaK_ETP(:,1) Notchfilter_lnCaK_ETP.data ./ std(Notchfilter_lnCaK_ETP.data)];

Interval_23_5 = [0.018 0.028];
Filter_obl_lnCaK_ETP_raw = idealfilter(Notchfilter_lnCaK_ETP,Interval_23_5,'pass');
Filter_obl_lnCaK_ETP = [lnCaK_ETP(:,1) Filter_obl_lnCaK_ETP_raw.data ./ std(Filter_obl_lnCaK_ETP_raw.data)];
Filter_obl_lnCaK_ETP_amp = [Filter_obl_lnCaK_ETP(:,1) abs(hilbert(Filter_obl_lnCaK_ETP(:,2)))];

Interval_9_3 = [0.006 0.012];
Filter_eccS_lnCaK_ETP_raw = idealfilter(Notchfilter_lnCaK_ETP,Interval_9_3,'pass');
Filter_eccS_lnCaK_ETP = [lnCaK_ETP(:,1) Filter_eccS_lnCaK_ETP_raw.data ./ std(Filter_eccS_lnCaK_ETP_raw.data)];
Filter_eccS_lnCaK_ETP_amp = [Filter_eccS_lnCaK_ETP(:,1) abs(hilbert(Filter_eccS_lnCaK_ETP(:,2)))];

Interval_25_05 = [0.002 0.003];
Filter_eccL_lnCaK_ETP_raw = idealfilter(Notchfilter_lnCaK_ETP,Interval_25_05,'pass');
Filter_eccL_lnCaK_ETP = [lnCaK_ETP(:,1) Filter_eccL_lnCaK_ETP_raw.data ./ std(Filter_eccL_lnCaK_ETP_raw.data)];

timeseries_lnCaK = timeseries(Filter_obl_lnCaK_ETP_amp(:,2));
Interval_obl_amp = [0.005 0.007];
lnCaK_obl_amp_filter_raw = idealfilter(timeseries_lnCaK,Interval_obl_amp,'pass');
lnCaK_obl_amp_filter = [Filter_obl_lnCaK_ETP_amp(:,1) lnCaK_obl_amp_filter_raw.data];

%% Load and Prepare d18O U1406
load('U1406_Isotopes.mat');
U1406_d18O_ETP_raw = [interp1(Tuning_ETP(:,1),Tuning_ETP(:,2),U1406_Isotopes(:,1),'linear','extrap') U1406_Isotopes(:,4)];
U1406_d18O_ETP_raw = U1406_d18O_ETP_raw(~isnan(U1406_d18O_ETP_raw(:,1)),:);
[~,a,~] = unique(U1406_d18O_ETP_raw(:,1));
U1406_d18O_ETP_raw = U1406_d18O_ETP_raw(a,:);

U1406_d18O_ETP = [21.772:0.001:26.349;interp1(U1406_d18O_ETP_raw(:,1),U1406_d18O_ETP_raw(:,2),21.772:0.001:26.349,'linear','extrap')]';
Data_notchfilter_d18O_ETP = timeseries(U1406_d18O_ETP(:,2));
Interval_notch_ETP = [0 0.001];
Notchfilter_d18O_ETP = idealfilter(Data_notchfilter_d18O_ETP,Interval_notch_ETP,'notch');
Notch_d18O_ETP = [U1406_d18O_ETP(:,1) Notchfilter_d18O_ETP.data];
Notch_d18O_ETP_plot = [U1406_d18O_ETP(:,1) Notchfilter_d18O_ETP.data ./ std(Notchfilter_d18O_ETP.data)];

Interval_23_5 = [0.018 0.028];
Filter_obl_U1406_d18O_raw = idealfilter(Notchfilter_d18O_ETP,Interval_23_5,'pass');
Filter_obl_U1406_d18O = [U1406_d18O_ETP(:,1) Filter_obl_U1406_d18O_raw.data ./ std(Filter_obl_U1406_d18O_raw.data)];
Filter_obl_U1406_d18O_amp = [Filter_obl_U1406_d18O(:,1) abs(hilbert(Filter_obl_U1406_d18O(:,2)))];

Interval_9_3 = [0.006 0.012];
Filter_eccS_U1406_d18O_raw = idealfilter(Notchfilter_d18O_ETP,Interval_9_3,'pass');
Filter_eccS_U1406_d18O = [U1406_d18O_ETP(:,1) Filter_eccS_U1406_d18O_raw.data ./ std(Filter_eccS_U1406_d18O_raw.data)];
Filter_eccS_U1406_d18O_amp = [Filter_eccS_U1406_d18O(:,1) abs(hilbert(Filter_eccS_U1406_d18O(:,2)))];

Interval_25_05 = [0.002 0.003];
Filter_eccL_U1406_d18O_raw = idealfilter(Notchfilter_d18O_ETP,Interval_25_05,'pass');
Filter_eccL_U1406_d18O = [U1406_d18O_ETP(:,1) Filter_eccL_U1406_d18O_raw.data ./ std(Filter_eccL_U1406_d18O_raw.data)];

timeseries_d18O = timeseries(Filter_obl_U1406_d18O_amp(:,2));
Interval_obl_amp = [0.005 0.007];
d18O_obl_amp_filter_raw = idealfilter(timeseries_d18O,Interval_obl_amp,'pass');
d18O_obl_amp_filter = [Filter_obl_U1406_d18O_amp(:,1) d18O_obl_amp_filter_raw.data];

%% Load La2004 Data
load('La2004_raw.mat');

La2004_age = La2004_ecc(:,1);

Data_obl = timeseries(La2004_obl(:,2) - mean(La2004_obl(:,2)));
Interval_23_5 = [0.018 0.028];
Filter_obl_23_5_raw = idealfilter(Data_obl,Interval_23_5,'pass');
Filter_obl_23_5 = [La2004_age(:,1) Filter_obl_23_5_raw.data ./ std(Filter_obl_23_5_raw.data)];
Filter_obl_23_5_amp = [La2004_age(:,1) abs(hilbert(Filter_obl_23_5(:,2)))];

Data_ecc = timeseries(La2004_ecc(:,2));
Interval_9_3 = [0.006 0.012];
Filter_eccS_raw = idealfilter(Data_ecc,Interval_9_3,'pass');
Filter_eccS = [La2004_age(:,1) Filter_eccS_raw.data ./ std(Filter_eccS_raw.data)];
Filter_eccS_amp = [La2004_age(:,1) abs(hilbert(Filter_eccS(:,2)))];

Interval_25_05 = [0.002 0.003];
Filter_eccL_raw = idealfilter(Data_ecc,Interval_25_05,'pass');
Filter_eccL = [La2004_age(:,1) Filter_eccL_raw.data ./ std(Filter_eccL_raw.data)];

timeseries_La2004obl = timeseries(Filter_obl_23_5_amp(:,2));
Interval_obl_amp = [0.005 0.007];
La2004_obl_amp_filter_raw = idealfilter(timeseries_La2004obl,Interval_obl_amp,'pass');
La2004_obl_amp_filter = [Filter_obl_23_5_amp(:,1) La2004_obl_amp_filter_raw.data];

%% Wavelet Analysis CaK
d_CaK = Notch_lnCaK_ETP;
[d_CaK,dt_CaK]=formatts(d_CaK);

variance_CaK = std(d_CaK(:,2)).^2;
n_CaK=size(d_CaK,1);
sigma2_CaK=var(d_CaK(:,2));
varargin = {};

%----------default arguments for the wavelet transform-----------
Args_CaK=struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/12, ...    % TVP: Changed to 96!!!! Was 12 sub-octaves per octave 
    'S0',2*dt_CaK,...    % this says start at a scale of 2 years
    'J1',[],...
    'Mother','Morlet', ...
    'MaxScale',[],...   %a more simple way to specify J1
    'AR1','auto');
Args_CaK=parseArgs(varargin,Args_CaK,{'BlackandWhite'});
if isempty(Args_CaK.J1)
    if isempty(Args_CaK.MaxScale)
        Args_CaK.MaxScale=(n_CaK*.17)*2*dt_CaK; %automaxscale
    end
    Args_CaK.J1=round(log2(Args_CaK.MaxScale/Args_CaK.S0)/Args_CaK.Dj);
end

if strcmpi(Args_CaK.AR1,'auto')
    Args_CaK.AR1=ar1nv(d_CaK(:,2));
    
    if any(isnan(Args_CaK.AR1))
        error('Automatic AR1 estimation failed. Specify it manually (use arcov or arburg).')
    end
end

[wave_CaK,period_CaK,scale_CaK,coi_CaK] = wavelet(d_CaK(:,2),dt_CaK,Args_CaK.Pad,Args_CaK.Dj,Args_CaK.S0,Args_CaK.J1,Args_CaK.Mother);

t_CaK = round(d_CaK(:,1)*1000)/1000;
power_CaK = (abs(wave_CaK)).^2 ;        % compute wavelet power spectrum
signif_CaK = wave_signif(1.0,dt_CaK,scale_CaK,0,Args_CaK.AR1,-1,-1,Args_CaK.Mother);
sig95_CaK = (signif_CaK')*(ones(1,n_CaK));  % expand signif --> (J+1)x(N) array
sig95_CaK = power_CaK ./ (sigma2_CaK*sig95_CaK);
Yticks_CaK = 2.^(fix(log2(min(period_CaK))):fix(log2(max(period_CaK))));

% Global Wavelet Power Spectrum
global_ws_CaK = mean(power_CaK,2)';   % time-average over all times
dof_CaK = n_CaK; %- scale_CaK;  % the -scale corrects for padding at edges
global_signif_CaK_95 = wave_signif(variance_CaK,dt_CaK,scale_CaK,1,Args_CaK.AR1,0.95,dof_CaK,Args_CaK.Mother);
global_signif_CaK_80 = wave_signif(variance_CaK,dt_CaK,scale_CaK,1,Args_CaK.AR1,0.80,dof_CaK,Args_CaK.Mother);

% Bias rectification
powers_CaK = zeros(size(power_CaK));
for k = 1:length(scale_CaK)
    powers_CaK(k,:) = power_CaK(k,:)/scale_CaK(k);
end
global_wss_CaK = global_ws_CaK ./ scale_CaK;
global_signifs_CaK_95 = global_signif_CaK_95 ./ scale_CaK;
global_signifs_CaK_80 = global_signif_CaK_80 ./ scale_CaK;

%% Wavelet Analysis d18O
d_d18O = Notch_d18O_ETP;
[d_d18O,dt_d18O]=formatts(d_d18O);

variance_d18O = std(d_d18O(:,2)).^2;
n_d18O = size(d_d18O,1);
sigma2_d18O = var(d_d18O(:,2));
varargin = {};

%----------default arguments for the wavelet transform-----------
Args_d18O = struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/12, ...    % TVP: Changed to 96!!!! Was 12 sub-octaves per octave 
    'S0',2*dt_d18O,...    % this says start at a scale of 2 years
    'J1',[],...
    'Mother','Morlet', ...
    'MaxScale',[],...   %a more simple way to specify J1
    'AR1','auto');
Args_d18O = parseArgs(varargin,Args_d18O,{'BlackandWhite'});
if isempty(Args_d18O.J1)
    if isempty(Args_d18O.MaxScale)
        Args_d18O.MaxScale = (n_d18O*.17)*2*dt_d18O; %automaxscale
    end
    Args_d18O.J1 = round(log2(Args_d18O.MaxScale/Args_d18O.S0)/Args_d18O.Dj);
end

if strcmpi(Args_d18O.AR1,'auto')
    Args_d18O.AR1 = ar1nv(d_d18O(:,2));
    
    if any(isnan(Args_d18O.AR1))
        error('Automatic AR1 estimation failed. Specify it manually (use arcov or arburg).')
    end
end

[wave_d18O,period_d18O,scale_d18O,coi_d18O] = wavelet(d_d18O(:,2),dt_d18O,Args_d18O.Pad,Args_d18O.Dj,Args_d18O.S0,Args_d18O.J1,Args_d18O.Mother);

t_d18O = round(d_d18O(:,1)*1000)/1000;
power_d18O = (abs(wave_d18O)).^2 ;        % compute wavelet power spectrum
signif_d18O = wave_signif(1.0,dt_d18O,scale_d18O,0,Args_d18O.AR1,-1,-1,Args_d18O.Mother);
sig95_d18O = (signif_d18O')*(ones(1,n_d18O));  % expand signif --> (J+1)x(N) array
sig95_d18O = power_d18O ./ (sigma2_d18O*sig95_d18O);
Yticks_d18O = 2.^(fix(log2(min(period_d18O))):fix(log2(max(period_d18O))));

% Global Wavelet Power Spectrum
global_ws_d18O = mean(power_d18O,2)';   % time-average over all times
dof_d18O = n_d18O - scale_d18O;  % the -scale corrects for padding at edges
global_signif_d18O_95 = wave_signif(variance_d18O,dt_d18O,scale_d18O,1,Args_d18O.AR1,0.95,dof_d18O,Args_d18O.Mother);
global_signif_d18O_80 = wave_signif(variance_d18O,dt_d18O,scale_d18O,1,Args_d18O.AR1,0.80,dof_d18O,Args_d18O.Mother);

% Bias rectification
powers_d18O = zeros(size(power_d18O));
for k = 1:length(scale_d18O)
    powers_d18O(k,:) = power_d18O(k,:)/scale_d18O(k);
end
global_wss_d18O = global_ws_d18O ./ scale_d18O;
global_signifs_d18O_95 = global_signif_d18O_95 ./ scale_d18O;
global_signifs_d18O_80 = global_signif_d18O_80 ./ scale_d18O;

%% Main figure
figure1 = figure('RendererMode','manual',...
    'Renderer','Painters',...
    'InvertHardCopy','off');

colormap('jet');

dcmObj = datacursormode;
set(dcmObj,'UpdateFcn',@custom_cursor);

% Data ln(Ca/K)
fig(01) = axes('Parent',figure1,...
    'position',pos01);
hold(fig(01),'all');
p1(01) = plot(Notch_lnCaK_ETP_plot(:,1),Notch_lnCaK_ETP_plot(:,2));

% Data ln(Ca/K)
fig(02) = axes('Parent',figure1,...
    'position',pos02);
hold(fig(02),'all');
p2(01) = plot(Notch_d18O_ETP_plot(:,1),-Notch_d18O_ETP_plot(:,2));

% Obliquity filters
fig(03) = axes('Parent',figure1,...
    'position',pos03);
hold(fig(03),'all');
p3(01) = plot(Filter_obl_23_5(:,1),Filter_obl_23_5(:,2));
p3(02) = plot(Filter_obl_lnCaK_ETP(:,1),Filter_obl_lnCaK_ETP(:,2));

% Obliquity amplitude
fig(04) = axes('Parent',figure1,...
    'position',pos04);
hold(fig(04),'all');
p4(01) = plot(La2004_obl_amp_filter(:,1),La2004_obl_amp_filter(:,2));
p4(02) = plot(lnCaK_obl_amp_filter(:,1),lnCaK_obl_amp_filter(:,2));

% short eccentricity filters
fig(05) = axes('Parent',figure1,...
    'position',pos05);
hold(fig(05),'all');
p5(01) = plot(Filter_eccS(:,1),Filter_eccS(:,2));
p5(02) = plot(Filter_eccS_U1406_d18O(:,1),-Filter_eccS_U1406_d18O(:,2));

% long eccentricity filters
fig(6) = axes('Parent',figure1,...
    'position',pos06);
hold(fig(6),'all');
p6(01) = plot(Filter_eccL(:,1),Filter_eccL(:,2));
p6(02) = plot(Filter_eccL_U1406_d18O(:,1),-Filter_eccL_U1406_d18O(:,2));

% Wavelet ln(Ca/K)
fig(7) = axes('Parent',figure1,...
    'position',pos07);
hold(fig(7),'all');
H=imagesc(t_CaK,log2(period_CaK),(log2(abs(powers_CaK/sigma2_CaK))));
colour_wavelet = reshape(log2(abs(powers_CaK/sigma2_CaK)),1,[]);
set(fig(7),'clim',[-1 11]);

[c,h] = contour(t_CaK,log2(period_CaK),(sig95_CaK),[1 1],'k');
pw(1) = plot(t_CaK,log2(coi_CaK));
for k = round(21:0.040:27,3)
    line([k k],[log2(coi_CaK(t_CaK == k)) 1],...
        'color','k','linewidth',0.5);
end

% Subplot Global Power Spectrum CaK
fig(8) = axes('Parent',figure1,...
    'Position',pos08);
hold(fig(8),'all');
p8(01) = plot(global_wss_CaK,log2(period_CaK));
p8(02) = plot(global_signifs_CaK_80,log2(period_CaK));
p8(03) = plot(global_signifs_CaK_95,log2(period_CaK));

% Wavelet d18O
fig(9) = axes('Parent',figure1,...
    'position',pos09);
hold(fig(9),'all');
H=imagesc(t_d18O,log2(period_d18O),(log2(abs(powers_d18O/sigma2_d18O))));
colour_wavelet = reshape(log2(abs(powers_d18O/sigma2_d18O)),1,[]);
set(fig(9),'clim',[1 10]);

[c,h] = contour(t_d18O,log2(period_d18O),(sig95_d18O),[1 1],'k');
pw(2) = plot(t_d18O,log2(coi_d18O));
for k = round(21:0.040:27,3)
    line([k k],[log2(coi_d18O(t_d18O == k)) 1],...
        'color','k','linewidth',0.5);
end

% Subplot Global Power Spectrum d18O
fig(10) = axes('Parent',figure1,...
    'Position',pos10);
hold(fig(10),'all');
p10(01) = plot(global_wss_d18O,log2(period_d18O));
p10(02) = plot(global_signifs_d18O_80,log2(period_d18O));
p10(03) = plot(global_signifs_d18O_95,log2(period_d18O));


%% Properties figures
set(fig(7:10),'box','on',...
    'xgrid','on',...
    'ygrid','on');

set(fig([1:7 9]),'xtick',xpos_age(2:end-1),...
    'xlim',[min(xpos_age) max(xpos_age)],...
    'xminortick','on',...
    'yminortick','on');
set(fig(2:10),'xticklabel',{});
set(fig(2:6),'xcolor','w');
set(fig([1:6 8 10]),'color','none');

set(fig(1:2),'ylim',[-3 3],...
    'ytick',-3:3);
set(fig(3),'ylim',[-4 4],...
    'ytick',-2:2:2);
set(fig(4),'ylim',[-1 1],...
    'ytick',-0.5:0.5:0.5);
set(fig(5),'ylim',[-4 4],...
    'ytick',-2:2:2);
set(fig(6),'ylim',[-2.1 2.1]);

set(fig(7:10),'ydir','reverse',...
    'ylim',[-7 0],...
    'ytick',log2([0.02 0.041 0.110 0.405]),...
    'layer','top');
set(fig([8 10]),'yticklabel',{});
set(fig(8),'xlim',[0 1.6e4],...
    'xtick',0:0.4e4:1.6e4);
set(fig([7 9]),'ytickLabel',{'20','41','110','405'});
set(fig(10),'xlim',[0 10],...
    'xtick',0:2.5:10);

%% Properties Plots
set(pw(1:2),'color','k');
set([p2(1) p5(2) p6(2) p10(1)],'color','b');
set([p1(1) p3(2) p4(2) p8(1)],'color',[0.6 0.6 0]);
set([p3(1) p4(1) p5(1) p6(1)],'color','r');
%set([p5(3) p6(3)],'color',[0 0.5 0],'linestyle','--');

set([p8(2) p10(2)],'linestyle',':','color','k');
set([p8(3) p10(3)],'linestyle','--','color','k');

%% Axes Labels
ylabel(fig(01),{'Processed','CaCO3'}); % ylabel(fig(01),{'a)','Detrend','Standard','-d18O','CaCO3'});
ylabel(fig(02),{'Processed','d18O'}); 
ylabel(fig(03),{'Obl','Filters'}); % ylabel(fig(02),{'b)','Filtered','Standard','-d18O','CaCO3, Obl','23 +/- 5','cycles/Myr'});
ylabel(fig(04),{'Obl AM','of filters'}); % ylabel(fig(03),{'c)','Amplitudes','Filt+Stand','d18O','CaCO3, Obl','23 +/- 5','cycles/Myr'});
ylabel(fig(05),{'Short Ecc','Filters'}); % ylabel(fig(04),{'d)','Filtered','Standard','-d18O','CaCO3, Ecc','9 +/- 3','cycles/Myr'});
ylabel(fig(06),{'Long Ecc','Filters'}); % ylabel(fig(05),{'e)','Amplitudes','Filt+Stand','d18O','CaCO3, Ecc','2.5 +/- 0.5','cycles/Myr'});
ylabel(fig(07),{'CaCO3','Period (kyr/cycle)'}); % ylabel(fig(06),{'f) & g) CaCO3','Period (kyr/cycle)'});
ylabel(fig(09),{'d18O','Period (kyr/cycle)'}); % ylabel(fig(08),{'h) & i) d18O','Period (kyr/cycle)'});

xlabel(fig(01),'Age (Ma) - Laskar et al. (2004) and this paper');
xlabel(fig(08),{'Mean','Power','\rightarrow'});

%legend(fig(03),{'La2004','d18O','CaCO3'});

%% Final layout
set(figure1,'units','centimeters',...
    'position',[0 0 21 20],...
    'outerposition',[0 0 21 20],...
    'paperunits','centimeters',...
    'resize','on',...
    'papersize',[21 20],...
    'paperposition',[0 0 21 20],...
    'paperpositionmode','manual');

cd(original_path);