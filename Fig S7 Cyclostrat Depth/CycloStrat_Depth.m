clear all
close all

original_path = cd();

ypos_depth = 35:5:175;
ypos_age = 21:0.2:27.2;

pos01 = [0.0700 0.015 0.020 0.855]; %Magstrat Depth
pos02 = [0.0950 0.015 0.060 0.855]; %Sed Rate Depth
pos03 = [0.1600 0.015 0.120 0.855]; %Raw ln(Ca/K) data
pos04 = [0.2850 0.015 0.210 0.855]; %Wavelet ln(Ca/K) Depth
pos05 = [pos04(1) 0.870 pos04(3) 0.09]; %Wavelet Global Power Spectrum ln(Ca/K)
pos06 = [0.5000 0.015 0.120 0.855]; %'Obl' Filtered ln(Ca/K) Depth
% Leave here a little space for correlation lines
pos07 = [0.6950 0.015 0.120 0.855]; %Obliquity La2004
pos08 = [0.8100 0.015 0.120 0.855]; %Eccentricity La2004

%% Load and Prepare Data ln(Ca/K)
lnCaK_data = xlsread('Data_lnCaK_depth.xlsx');

lnCaK_raw = [lnCaK_data(:,1) 5.7525 .* exp(0.6079 .* lnCaK_data(:,2))];

lnCaK_smooth = [lnCaK_raw(:,1) smooth(lnCaK_raw(:,2),'moving',5)];
lnCaK_interp = [34.66:0.02:174.72; interp1(lnCaK_raw(:,1),lnCaK_raw(:,2),34.66:0.02:174.72,'linear','extrap')]';

Data_notchfilter_lnCaK = timeseries(lnCaK_interp(:,2));
Interval_notch = [0 0.0008]; %Periods larger than 25 m -->>> 0.0004 cycle/cm = 0.0008 cycle/(2cm).
Notchfilter = idealfilter(Data_notchfilter_lnCaK,Interval_notch,'notch');
lnCaK_detrend = [lnCaK_interp(:,1) Notchfilter.data];

lnCaK_Int1 = lnCaK_detrend(1:3268,:); %Int1 top-100 m
lnCaK_Int2 = lnCaK_detrend(2518:5268,:); %Int2 85-140 m
lnCaK_Int3 = lnCaK_detrend(5018:end,:); %Int3 135-end m

Data_Filter_Int1 = timeseries(lnCaK_Int1(:,2));
Filter_obl_Int1_range = [0.012 0.0192]; %Depth scale is every 2 cm!! 0.78 +/- 0.18 cycles/m is 0.012 0.0192 cycles/(2cm)!!
Filter_obl_Int1_raw = idealfilter(Data_Filter_Int1,Filter_obl_Int1_range,'pass');
Filter_obl_Int1 = [lnCaK_Int1(:,1) Filter_obl_Int1_raw.data];
Filter_ecc_Int1_range = [0.0048 0.0092]; %Depth scale is every 2 cm!! 0.35 +/- 0.11 cycles/m is 0.0048 0.0092 cycles/(2cm)!!
Filter_ecc_Int1_raw = idealfilter(Data_Filter_Int1,Filter_ecc_Int1_range,'pass');
Filter_ecc_Int1 = [lnCaK_Int1(:,1) Filter_ecc_Int1_raw.data];
Filter_obl_Int1_amp = [Filter_obl_Int1(:,1) abs(hilbert(Filter_obl_Int1(:,2)))];
Filter_ecc_Int1_amp = [Filter_ecc_Int1(:,1) abs(hilbert(Filter_ecc_Int1(:,2)))];

Data_Filter_Int2 = timeseries(lnCaK_Int2(:,2));
Filter_obl_Int2_range = [0.0192 0.0308]; %Depth scale is every 2 cm!! 1.25 +/- 0.29 cycles/m is 0.0192 0.0308 cycles/(2cm)!!
Filter_obl_Int2_raw = idealfilter(Data_Filter_Int2,Filter_obl_Int2_range,'pass');
Filter_obl_Int2 = [lnCaK_Int2(:,1) Filter_obl_Int2_raw.data];
Filter_ecc_Int2_range = [0.0094 0.0146]; %Depth scale is every 2 cm!! 0.60 +/- 0.13 cycles/m is 0.0094 0.0146 cycles/(2cm)!!
Filter_ecc_Int2_raw = idealfilter(Data_Filter_Int2,Filter_ecc_Int2_range,'pass');
Filter_ecc_Int2 = [lnCaK_Int2(:,1) Filter_ecc_Int2_raw.data];
Filter_obl_Int2_amp = [Filter_obl_Int2(:,1) abs(hilbert(Filter_obl_Int2(:,2)))];
Filter_ecc_Int2_amp = [Filter_ecc_Int2(:,1) abs(hilbert(Filter_ecc_Int2(:,2)))];

Data_Filter_Int3 = timeseries(lnCaK_Int3(:,2));
Filter_obl_Int3_range = [0.0148 0.0204]; %Depth scale is every 2 cm!! 0.88 +/- 0.14 cycles/m is 0.0148 0.0204 cycles/(2cm)!!
Filter_obl_Int3_raw = idealfilter(Data_Filter_Int3,Filter_obl_Int3_range,'pass');
Filter_obl_Int3 = [lnCaK_Int3(:,1) Filter_obl_Int3_raw.data];
Filter_ecc_Int3_range = [0.0034 0.0098]; %Depth scale is every 2 cm!! 0.33 +/- 0.16 cycles/m is 0.0034 0.0098 cycles/(2cm)!!
Filter_ecc_Int3_raw = idealfilter(Data_Filter_Int3,Filter_ecc_Int3_range,'pass');
Filter_ecc_Int3 = [lnCaK_Int3(:,1) Filter_ecc_Int3_raw.data];
Filter_obl_Int3_amp = [Filter_obl_Int3(:,1) abs(hilbert(Filter_obl_Int3(:,2)))];
Filter_ecc_Int3_amp = [Filter_ecc_Int3(:,1) abs(hilbert(Filter_ecc_Int3(:,2)))];

%% Load Sed Rate U1406 Depth
[AgeModel_U1406_raw,AgeModel_U1406_text,~] = xlsread('AgeModel_U1406.xlsx');
AgeModel_U1406_isnan = ~isnan(AgeModel_U1406_raw(:,1));

Reversals = AgeModel_U1406_text(19:45,1);
AgeModel_U1406_depth = AgeModel_U1406_raw(AgeModel_U1406_isnan,1);
AgeModel_U1406_depth = AgeModel_U1406_depth(2:28);
Age_U1406_GTS2012 = AgeModel_U1406_raw(AgeModel_U1406_isnan,3);
Age_U1406_GTS2012 = Age_U1406_GTS2012(2:28);

SedRate_calc_GTS2012 = 1/10 * (AgeModel_U1406_depth(1:length(AgeModel_U1406_depth)-1,1) - AgeModel_U1406_depth(2:length(AgeModel_U1406_depth),1)) ./ (Age_U1406_GTS2012(1:length(Age_U1406_GTS2012)-1,1) - Age_U1406_GTS2012(2:length(Age_U1406_GTS2012),1));
SedRate_GTS2012_age_unsorted = [SedRate_calc_GTS2012 Age_U1406_GTS2012(1:length(Age_U1406_GTS2012)-1,1); SedRate_calc_GTS2012 Age_U1406_GTS2012(2:length(Age_U1406_GTS2012),1)-1e-6];
SedRate_GTS2012_depth_unsorted = [SedRate_calc_GTS2012 AgeModel_U1406_depth(1:length(AgeModel_U1406_depth)-1,1); SedRate_calc_GTS2012 AgeModel_U1406_depth(2:length(AgeModel_U1406_depth),1)-1e-6];
SedRate_GTS2012_age = sortrows(SedRate_GTS2012_age_unsorted,2);
SedRate_GTS2012_depth = sortrows(SedRate_GTS2012_depth_unsorted,2);

%% Load tuning spreadsheet
Tuning = xlsread('Tuning_TiePoints.xlsx');

%% Load U-channel Magstrat
load('Magstrat_U1406.mat');
U1406_Magstrat_cdata(U1406_Magstrat_cdata == 0.5) = 0.75;
U1406_Magstrat_cdata(U1406_Magstrat_cdata == 0) = 0.5;

%% Load La2004
load('La2004_raw.mat');
La2004_age = La2004_ecc(:,1);

env_obl_raw =  abs(hilbert(La2004_obl(:,2) - mean(La2004_obl(:,2))));
env_obl = [La2004_age env_obl_raw + mean(La2004_obl(:,2))];

Data_ecc = timeseries(La2004_ecc(:,2));
Interval_25_05 = [0.002 0.003];
Interval_04_02 = [0.0002 0.0006];
Filter_eccL400_raw = idealfilter(Data_ecc,Interval_25_05,'pass');
Filter_eccL24_raw = idealfilter(Data_ecc,Interval_04_02,'pass');
Filter_eccL400 = [La2004_age(:,1) Filter_eccL400_raw.data ./ std(Filter_eccL400_raw.data)];
Filter_eccL24 = [La2004_age(:,1) Filter_eccL24_raw.data ./ std(Filter_eccL24_raw.data)];


%% Wavelet Analysis CaK
d_CaK = lnCaK_detrend;
[d_CaK,dt_CaK]=formatts(d_CaK);

variance_CaK = std(d_CaK(:,2)).^2; 
n_CaK=size(d_CaK,1);
sigma2_CaK=var(d_CaK(:,2));
varargin = {};

%----------default arguments for the wavelet transform-----------
Args_CaK=struct('Pad',1,...      % pad the time series with zeroes (recommended)
    'Dj',1/48, ...    % TVP: Changed to 96!!!! Was 12 sub-octaves per octave 
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

t_CaK = round(d_CaK(:,1)*100)/100;
power_CaK = (abs(wave_CaK)).^2 ;        % compute wavelet power spectrum
signif_CaK = wave_signif(1.0,dt_CaK,scale_CaK,0,Args_CaK.AR1,-1,-1,Args_CaK.Mother);
sig95_CaK = (signif_CaK')*(ones(1,n_CaK));  % expand signif --> (J+1)x(N) array
sig95_CaK = power_CaK ./ (sigma2_CaK*sig95_CaK);

% Global Wavelet Power Spectrum
global_ws_CaK = mean(power_CaK,2)'; %variance_CaK .* (sum(power_CaK,2) ./ n_CaK);   % time-average over all times
dof_CaK = n_CaK - scale_CaK;  % the -scale corrects for padding at edges
[global_signif_CaK_95,fft_95] = wave_signif(variance_CaK,dt_CaK,scale_CaK,1,Args_CaK.AR1,0.95,dof_CaK,Args_CaK.Mother);
[global_signif_CaK_80,fft_80] = wave_signif(variance_CaK,dt_CaK,scale_CaK,1,Args_CaK.AR1,0.80,dof_CaK,Args_CaK.Mother);

[signif,fft] = wave_signif(variance_CaK,dt_CaK,scale_CaK,1,Args_CaK.AR1,0.95,dof_CaK,Args_CaK.Mother);


% Bias rectification
powers_CaK = zeros(size(power_CaK));
for k = 1:length(scale_CaK)
    powers_CaK(k,:) = power_CaK(k,:)/scale_CaK(k);
end
global_wss_CaK = global_ws_CaK ./ scale_CaK;
global_signifs_CaK_95 = global_signif_CaK_95 ./ scale_CaK;
global_signifs_CaK_80 = global_signif_CaK_80 ./ scale_CaK;


%% Generate figure
figure1 = figure('RendererMode','manual',...
    'Renderer','Painters',...
    'InvertHardCopy','off');
cmap1 = colormap('jet');
cmap2 = colormap('gray');
cmap = [cmap1;cmap2];
colormap(cmap);

dcmObj = datacursormode;
set(dcmObj,'UpdateFcn',@custom_cursor);

% Sed Rate Depth
fig(01) = axes('Parent',figure1,...
    'position',pos01);
hold(fig(01),'all');
pa(01) = patch('Vertices',U1406_Magstrat_depth,...
    'Faces',U1406_Magstrat_faces,...
    'FaceColor','flat',...
    'FaceVertexCData',U1406_Magstrat_cdata,...
    'CDataMapping','scaled',...
    'edgecolor','none');
pa(02) = patch('vertices',[0 min(ypos_depth);...
    1 min(ypos_depth);...
    1 max(ypos_depth);...
    0 max(ypos_depth)],...
    'faces',[1 2 3 4],...
    'facecolor','none',...
    'cdatamapping','scaled',...
    'edgecolor','k');

te(1) = text(0.5,43,'C6AAr.1r');
te(2) = text(0.5,76,'C6Br');
te(3) = text(0.5,102,'C6Cr');
te(4) = text(0.5,118,'C7r');
te(5) = text(0.5,152,'C8r');

te(6) = text(0.5,67,'C6Bn.2n');
te(7) = text(0.5,138,'C8n');
te(8) = text(0.5,166,'C9n');

% Revised Magstrat (van Peer et al., 2017)
fig(02) = axes('parent',figure1,...
    'position',pos02);
hold(fig(02),'all');
pl(1) = plot(SedRate_GTS2012_depth(:,1),SedRate_GTS2012_depth(:,2));

% Subplot 3 XRF Data Raw
fig(3) = axes('Parent',figure1,...
    'Position',pos03);
hold(fig(3),'all');
pl(2) = plot(lnCaK_raw(:,2),lnCaK_raw(:,1));
e(1) = errorbar(28,67.5,32.5);
e(2) = errorbar(23,112.5,27.5);
e(3) = errorbar(58,154.86,19.86);

% Wavelet XRF Data
fig(4) = axes('Parent',figure1,...
    'Position',pos04);
hold(fig(4),'all');
H=imagesc(log2(period_CaK),t_CaK,(log2(abs(powers_CaK/sigma2_CaK)))');
colour_wavelet = reshape(log2(abs(powers_CaK/sigma2_CaK)),1,[]);

clim_new = [-6 21];
set(fig(04),'clim',clim_new);

[c,h] = contour(log2(period_CaK),t_CaK,(sig95_CaK)',[1 1],'k');
%set(h,'linewidth',2);
tt=[t_CaK([1 1])-dt_CaK*.5;t_CaK;t_CaK([end end])+dt_CaK*.5];
pw(1) = plot(log2(coi_CaK),t_CaK);
for k = 36:1:174
    line([log2(coi_CaK(t_CaK == k)) 100],[k k],...
        'color','k','linewidth',0.5);
end

% Subplot 5 Global Power Spectrum CaK
fig(5) = axes('Parent',figure1,...
    'Position',pos05);
hold(fig(5),'all');
pw(2) = plot(log2(period_CaK),global_wss_CaK);
pw(3) = plot(log2(period_CaK),global_signifs_CaK_80);
pw(4) = plot(log2(period_CaK),global_signifs_CaK_95);


% Subplot 6 XRF Data Detrended and Filtered
fig(6) = axes('Parent',figure1,...
    'Position',pos06);
hold(fig(6),'all');
pl(05) = plot(lnCaK_detrend(:,2),lnCaK_detrend(:,1));
pl(06) = plot(2*Filter_obl_Int1(:,2),Filter_obl_Int1(:,1));
pl(07) = plot(2*Filter_obl_Int1_amp(:,2),Filter_obl_Int1_amp(:,1));
pl(08) = plot(2*Filter_obl_Int2(:,2),Filter_obl_Int2(:,1));
pl(09) = plot(2*Filter_obl_Int2_amp(:,2),Filter_obl_Int2_amp(:,1));
pl(10) = plot(2*Filter_obl_Int3(:,2),Filter_obl_Int3(:,1));
pl(11) = plot(2*Filter_obl_Int3_amp(:,2),Filter_obl_Int3_amp(:,1));

% Subplot 7 La2004 Obliquity
fig(07) = axes('Parent',figure1,...
    'Position',pos07);
hold(fig(07),'all');
pl(12) = plot(La2004_obl(:,2),La2004_age);
pl(13) = plot(env_obl(:,2),La2004_age);

% Subplot 8 La2004 Eccentricity
fig(08) = axes('Parent',figure1,...
    'Position',pos08);
hold(fig(08),'all');
pl(14) = plot(100*La2004_ecc(:,2),La2004_age);
pl(15) = plot(Filter_eccL400(:,2)+4.5,Filter_eccL400(:,1));
pl(16) = plot(Filter_eccL24(:,2)+5.5,Filter_eccL24(:,1));


%% Annotations of correlation lines
for k = 1:length(Tuning)
    an = Tuning(k,1);
    line([20 25],[an an],...
        'parent',fig(06),...
        'color','k');
    an2 = Tuning(k,2);
    line([22 22.5],[an2 an2],...
        'parent',fig(07),...
        'color','k');
    annotation(figure1,'line',[pos06(1) + pos06(3) pos07(1)],...
        [(1-(an-min(ypos_depth)) ./ (max(ypos_depth) - min(ypos_depth))) ...
        .* pos06(4)  + pos06(2)...
        (1-(an2-min(ypos_age)) ./ (max(ypos_age) - min(ypos_age))) ...
        .* pos07(4)  + pos07(2)],...
        'color','k');
    % Diagonal lines drawn between right side of plot 6 (0.62) and left 
    % side of plot 7 (0.695). First subtract 35 (top of plot in depth) and
    % 20.8 (top of plot in age), then devide by total length (140 m in
    % depth, 6.4 Myr in age). But y-axis is reversed, so 1-result is what
    % we need for plotting. Then multiply by total y-length of the plot and
    % move according to offset relative to the base of the page
end

%% Layout axes
for k = 2:5
    box(fig(k),'on');
end

set(fig([1:4 6]),...
    'ydir','reverse',...
    'ylim',[min(ypos_depth) max(ypos_depth)],...
    'ytick',ypos_depth);
set(fig(7:8),...
    'ydir','reverse',...
    'ylim',[min(ypos_age) max(ypos_age)],...
    'ytick',ypos_age);
set(fig(2:7),'yticklabel',{});
set(fig(8),'yaxislocation','right');

set(fig(7),'ycolor','w');
set(fig([2:3 5:8]),'color','none');

set(fig(2:5),'xgrid','on');
set(fig([1:3 6:8]),...
    'ygrid','on',...
    'yminortick','on');
set(fig([1 3 6:8]),'xminortick','on');

set(fig(01),'xtick',[],...
    'tickdir','out',...
    'ticklength',[0.004 0.004],...
    'clim',[0 1],...
    'xaxislocation','top');
xlabel(fig(01),'a)');
ylabel(fig(01),'Depth revised CCSF-A (m) - van Peer et al. (2017a)');

set(fig(02),'xlim',[0.5 4],...
    'xticklabel',{1 2 3 ''},...
    'xaxislocation','top');
xlabel(fig(02),{'b)','LSR','(cm/kyr)'});

set(fig(03),'xlim',[15 65],...
    'xtick',[20 40 60],...
    'xaxislocation','top');
xlabel(fig(03),{'c)','CaCO3'});

set(fig(04),'xlim',[log2(0.15) log2(40)],...
    'xtick',log2([0.25 0.5 1 2 4 8 16 32]),...
    'xticklabel',{},...
    'layer','top',...
    'ygrid','on');

set(fig(05),'xlim',get(fig(04),'xlim'),...
    'xtick',get(fig(04),'xtick'),...
    'xticklabel',{0.25 '' 1 2 4 8 16 32},...
    'xaxislocation','top',...
    'ylim',[0 600],...
    'ytick',0:150:600,...
    'yaxislocation','right',...
    'ygrid','on');
ylabel(fig(05),{'Mean Power \rightarrow'});
xlabel(fig(05),{'d) Period (m/cycle)'});

set(fig(6),'xlim',[-25 25],...
    'xtick',-20:10:20,...
    'xticklabel',{'' -10:10:20},...
    'xaxislocation','top');
xlabel(fig(06),{'f)','Filtered','CaCO3'});

set(fig(7),'xlim',[22 24.5],...
    'xtick',22:24,...
    'xticklabel',{22 23 24},...
    'xaxislocation','top');
xlabel(fig(7),{'h)','La2004','Obliquity','& envelope','(degrees)'});
set(fig(8),'xlim',[0 8],...
    'xtick',0:3:6,...
    'xaxislocation','top');
xlabel(fig(8),{'i)','La2004','Eccentricity','& filters','(%)'});
ylabel(fig(8),'Age (Ma) - Laskar et al. (2004)');

%% Layout plots
set([pw(1:2) pl([1 13 15])],'color','k');
set(pw(3),'linestyle',':','color','k');
set(pw(4),'linestyle','--','color','k');

set(pl([5 12 14]),'color',[0.3 0.3 0.3]);
set(pl(2),'color',[0.6 0.6 0]);
%set(pl([2 19]),'color',[0.5 0.5 0.5]);
set(pl([6 7]),'color','b');
set(pl([8 9]),'color','r');
set(pl([10 11]),'color',[0.0 0.8 0.8]);
set(e(1),'color','b');
set(e(2),'color','r');
set(e(3),'color',[0.0 0.8 0.8]);

%% Final layout
set(figure1,'units','centimeters',...
    'position',[0 0 21 29],...
    'outerposition',[0 0 21 29],...
    'paperunits','centimeters',...
    'resize','off',...
    'papersize',[21 29],...
    'paperposition',[0 0 21 29],...
    'paperpositionmode','manual');

cd(original_path);