% Amirreza Bahramani
% Advanced Neuroscience Course
% Sharif University
% June 2023
% Final Project
% Simulating Results of
% Fayaz, S., Fakharian, M. A. & Ghazizadeh, A. Stimulus presentation can...
% enhance spiking irregularity across subcortical and cortical regions.
% PLoS Comput Biol 18, e1010256 (2022).
% And

% close all;
clear;
clc;

addpath('Data_Analyzed\')
addpath(genpath('C:\Users\rimaz\Documents\Personal Documents\M.Sc. E.E. at Sharif\Courses\Semester 2\Advanced Neuroscience\Final Project\Code_and_Data\lib'))

%%
clc;
clear;

% Two SNR Samples from Good and Bad sessions
% load("Data_Raw\Data_Trans\myHVC.mat")
% load("Data_Raw\Data_Trans\stim_all.mat")
% load Data_Raw/SNR_Sample_Bad.mat
% load Data_Raw/SNR_Sample_Good.mat
% load("Data_Raw\Data_Trans\MGd.mat")
% load("Data_Raw\Data_Trans\stim_all.mat")
load("Data_Raw\Data_Trans\IpN1.mat")

neuronID = 1;
raster_{1} = raster{neuronID};

% raster_ = raster;

% start/end times
params.t_start = -500;
params.t_end = 1500;
% sliding window and binwidth in ms
params.slide = 50;
params.bin_width = 300;

% CV2_Local
cv = var_cv(raster_,params);
params.kass_params.cv2_all_neu = cv.cv2_all_neu;
params.kass_params.K = 14;

% Var Decompostion
result = var_decom(raster_,params);
result.cv = cv;

% VarCE
result.varce = var_varce(raster_{1},params);

t_start = params.t_start;
t_end = params.t_end;
bin_width = params.bin_width;

raster_local = raster_{1};
n = numel(raster_local);

bin_size = 20;

t_rate = t_start:t_end-bin_size;

fr = zeros(n,length(t_rate));

spike_mat = cell2mat(cellfun(@(x) histcounts(x,t_start:t_end),...
    raster_local,'UniformOutput',false));

for i = 1:n
    perccount(i,numel(raster))
    fr_temp = spike_mat(i,:);
    fr(i,:) = movmean(fr_temp,bin_size,'Endpoints','discard')*1e3;
end

rate = mean(fr);
rate_sd = std(fr);
t_rate = t_start:t_end-bin_size;

result.raster = raster_;
result.rate.rate = rate;
result.rate.t_rate = t_rate;
result.rate.rate_trial = fr;
result.rate.rate_sd = rate_sd;

%% plot

T = [-150,300];

red = [255, 35, 71]/255;
blue = [0, 145, 207]/255;
orange = [255, 111, 58]/255;
skyblue = [0, 188, 236]/255;
green = [0, 148, 50]/255;
lightgreen = [163, 203, 56]/255;
viol = [142, 68, 173]/255;
pink = 	[217, 128, 250]/255;

k = 20;

raster = result.raster{1};
n = numel(raster);

if n<k
    k = n;
end

rng('default');
tr_inds = randsample(n,k);

r = 4;
c = 4;


%%
figure
set(gcf, 'Position', [13 360 1900 620])
% FR
subplot(2,4,[1,2]);
plotErrorbar(result.rate.t_rate,result.rate.rate,...
    result.rate.rate_sd,'b');
xline(0,'r--')
xlim([t_start+bin_width/2,t_end-bin_width/2]);
xticklabels({})
ylabel('Firing Rate [Hz]', 'Interpreter','latex')

subplot(2,4,[5,6]);
raster_ = raster(tr_inds);

plot_raster(raster_,'b',.5);
axis tight;
box off;
xline(0,'r--')
xlim([t_start+bin_width/2,t_end-bin_width/2]);
ylabel('Trial Num', 'Interpreter','latex')
xlabel('Time [ms]', 'Interpreter','latex')

% plot nRV
subplot(2,4,3)
set(gca,'colororder',[green;red])
yyaxis left
plotErrorbar(result.times,result.nRV.kass,result.nRV.kass_SE,green);
% plotErrorbar(result.times,result.VEC.kass,result.VEC.kass_SE,green);
ylabel('nRV', 'Interpreter','latex')

% plot FF
yyaxis right
plotErrorbar(result.times,result.FF,result.FF_SE,red);
ylabel('Fano Factor', 'Interpreter','latex')

xline(0,'r--')
xlim([t_start+bin_width/2,t_end-bin_width/2]);
xticklabels({})

% plot nSI
subplot(2,4,7)
set(gca,'colororder',[viol;lightgreen])
yyaxis left
plotErrorbar(result.times,result.nSI.kass,result.nSI.kass_SE,viol);
% plotErrorbar(times,res.(name).EVC.kass,res.(name).EVC.kass_SE,viol);
ylabel('$n\Psi$', 'Interpreter','latex')

% plot CV2
yyaxis right
plotErrorbar(result.times,...
    result.cv.cv2,result.cv.cv2_SE,lightgreen);
ylabel('Local CV$\textsuperscript{2}$', 'Interpreter','latex')

xline(0,'r--')
ylim([0,.8])
xlabel('Time [ms]', 'Interpreter','latex')
xlim([t_start+bin_width/2,t_end-bin_width/2]);


subplot(2,4,4)
% VEC VarCE
h1=plotErrorbar(result.varce.times,...
    double(result.varce.VEC.VarCE),double(result.varce.VEC.VarCE_SE),red);
hold on
% VEC Vinci
h2=plotErrorbar(result.times,result.VEC.kass,result.VEC.kass_SE,blue);
ylabel('VEC', 'Interpreter','latex')
xlim([t_start+bin_width/2,t_end-bin_width/2]);
legend ([h1, h2],'VarCE','Kass','Location','best', 'Interpreter','latex');
legend('boxoff');
xline(0,'r--','HandleVisibility','off')
xticklabels({})

subplot(2,4,8)
% EVC VarCE
h1=plotErrorbar(result.varce.times,...
    result.varce.nSI.VarCE .*result.FR*bin_width*1e-3,...
    result.varce.nSI.VarCE .*result.FR_SE*bin_width*1e-3,red);

hold on
% EVC Vinci
h2=plotErrorbar(result.times,result.EVC.kass,result.EVC.kass_SE,blue);
ylabel('EVC', 'Interpreter','latex')
xlim([t_start+bin_width/2,t_end-bin_width/2]); 
xline(0,'r--')
legend ([h1, h2],'VarCE','Kass','Location','best', 'Interpreter','latex');legend('boxoff');

sgtitle(['Single IpN Neuron'], 'Interpreter','latex')
