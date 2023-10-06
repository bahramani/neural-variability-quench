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

close all;
clear;
clc;

addpath('Data_Analyzed\')
addpath(genpath('C:\Users\rimaz\Documents\Personal Documents\M.Sc. E.E. at Sharif\Courses\Semester 2\Advanced Neuroscience\Final Project\Code_and_Data\lib'))

%% UPLOAD DATA
areaParam = struct();

% SNR
areaParam(1).name = 'SNR';
areaParam(1).T = [-100,50];  % not important 
areaParam(1).llim = 1.3;
areaParam(1).dlim = 0;
areaParam(1).align = 0;
areaParam(1).bin_width = 200; % important
areaParam(1).cv_lim = [.27,.47];

% vlPFC
areaParam(2).name = 'vlPFC';
areaParam(2).T = [-100,100];
areaParam(2).llim = 3;
areaParam(2).dlim = 0.5;
areaParam(2).align = 0;
areaParam(2).bin_width = 200;
areaParam(2).cv_lim = [.4,.8];

% IpN1
areaParam(3).name = 'IpN1';
areaParam(3).T = [-200,250];
areaParam(3).llim = 1.5;
areaParam(3).dlim = 0;
areaParam(3).align = 0;
areaParam(3).bin_width = 300;
areaParam(3).cv_lim = [0.2 0.4];

% IpN2
areaParam(4).name = 'IpN2';
areaParam(4).T = [-100,350];
areaParam(4).llim = 2;
areaParam(4).dlim = 0;
areaParam(4).align = 0;
areaParam(4).bin_width = 300;
areaParam(4).cv_lim = [0.2 0.4];

% L/CLM All
areaParam(5).name = 'Stim_10';
areaParam(5).T = [-200,500];
areaParam(5).llim = 5;
areaParam(5).dlim = 1;
areaParam(5).align = 0;
areaParam(5).bin_width = 300;
areaParam(5).cv_lim = [0.7,1.1];

% % L/CLM Stim
% areaParam(5).name = 'Stim_10';
% areaParam(5).T = [-200,500];
% areaParam(5).llim = 3;
% areaParam(5).dlim = 0.5;
% areaParam(5).align = 0;
% areaParam(5).bin_width = 50;
% areaParam(5).cv_lim = [.4,.8];

% MGd
areaParam(6).name = 'MGd';
areaParam(6).T = [-200,300];
areaParam(6).llim = 3;
areaParam(6).dlim = 0;
areaParam(6).align = 0;
areaParam(6).bin_width = 500;
areaParam(6).cv_lim = [.4,.8];

% VISam
areaParam(7).name = 'VISam';
areaParam(7).T = [-200,300];
areaParam(7).llim = 3;
areaParam(7).dlim = 0;
areaParam(7).align = 0;
areaParam(7).bin_width = 500;
areaParam(7).cv_lim = [0,1];

% VISp
areaParam(8).name = 'VISp';
areaParam(8).T = [-200,300];
areaParam(8).llim = 3;
areaParam(8).dlim = 0;
areaParam(8).align = 0;
areaParam(8).bin_width = 500;
areaParam(8).cv_lim = [0,1.5];

areaNames = {'SNR', 'vlPFC', 'V1', 'V2', 'MT', 'PMd', 'NCM'};

%% PARAMETERS
% figure
% set(gcf, 'WindowState', 'maximized')
% figName = 'Distribution of Different Biases';
% title(figName)
% if saveFig
%     saveas(gcf,['pics' figName '.png']);
% end
saveFig = false;
res = struct();

%%
red = [255, 35, 71]/255;
blue = [0, 159, 255]/255;
orange = [199, 96, 0]/255;
skyblue = [0, 188, 236]/255;
green = [0, 89, 3]/255;
lightgreen = [163, 203, 56]/255;
viol = [115, 0, 150]/255;
pink = 	[217, 128, 250]/255;
gray = 	.6*[1, 1, 1];

%% Population Response
close all;

for counterName = 5:5 %numel(areaNames)

    % parameters
    name = areaParam(counterName).name;
    T = areaParam(counterName).T;
    llim = areaParam(counterName).llim;
    dlim = areaParam(counterName).dlim;
    align = areaParam(counterName).align;
    bin_width = areaParam(counterName).bin_width;
    
    % load data
    load([name, '_Analyzed.mat']);
    res.(name) = results;
    clear results;

    numNeurons = length(res.(name).valid_index);
    times = res.(name).times - align;
    figure
    set(gcf, 'Position', [13 360 1900 620])
    
    % plot Rate at 10 ms bin
    t_FR = res.(name).rate.t - align;
    subplot(2,3,1)
    plot(t_FR,res.(name).rate.fr,'black','LineWidth',1.5)
    xline(0,'r--')
    % plot Rate at 300 ms bin (like FF)
    hold on;
    plotErrorbar(times,res.(name).FR,res.(name).FR_SE,gray);
    xlim([times(1),times(end)])
    xticklabels({})
    ylabel('Firing Rate [Hz]', 'Interpreter','latex')
    title(['Firing Rate of Population (n=', num2str(numNeurons), ')'], 'Interpreter','latex')
    
    % effect of time-bin          
    subplot(2,3,[3,6])  
    rectangle('Position',[.7*bin_width,-10,.3*bin_width,20],...
        'FaceColor',.9*[1,1,1]);  
    hl1 = plotErrorbar(res.(name).FFvT.T, res.(name).FFvT.pre.FF_mean,...
        res.(name).FFvT.pre.FF_SE,blue);
    hold on;
    hl2 = plotErrorbar(res.(name).FFvT.T, res.(name).FFvT.post.FF_mean,...
        res.(name).FFvT.post.FF_SE,orange);
    x = res.(name).FFvT.T; 
    plot(x,res.(name).FFvT.inter_pre+x/bin_width*res.(name).FFvT.sl_pre,'--',...
        'Color',blue);
    plot(x,res.(name).FFvT.inter_post+x/bin_width*res.(name).FFvT.sl_post,...
        '--','Color',orange);
    scatter(0,res.(name).FFvT.inter_pre,'filled','CData',blue);
    scatter(0,res.(name).FFvT.inter_post,'filled','CData',orange);
    errorbar(0,res.(name).FFvT.inter_pre,res.(name).FFvT.inter_pre_SE,...
        'Color',blue);
    errorbar(0,res.(name).FFvT.inter_post,res.(name).FFvT.inter_post_SE,...
        'Color',orange); 
    ylim([dlim,llim])
    xlabel('time-bin [ms]', 'Interpreter','latex')
    Ta = T - align;
    legend([hl1, hl2],{['Pre Stim. (',num2str(Ta(1)),'ms)'],...
        ['Post Stim. (',num2str(Ta(2)),'ms)']},'Box','off',...
        'Location','southeast', 'Interpreter','latex');


    % plot nRV
    subplot(2,3,2)
    set(gca,'colororder',[green;red])
    yyaxis left
    plotErrorbar(times,res.(name).nRV.kass,res.(name).nRV.kass_SE,green);
    % plotErrorbar(times,res.(name).VEC.kass,res.(name).VEC.kass_SE,green);
    ylabel('nRV', 'Interpreter','latex')

    % plot FF
    yyaxis right
    plotErrorbar(times,res.(name).FF,res.(name).FF_SE,red);
    ylabel('Fano Factor', 'Interpreter','latex')

    xline(0,'r--')
    xlim([times(1),times(end)])
    xticklabels({})
     
 
    % plot nRV
    subplot(2,3,4)
    set(gca,'colororder',[green;viol])
    yyaxis left
    plotErrorbar(times,res.(name).nRV.kass,res.(name).nRV.kass_SE,green);
    % plotErrorbar(times,res.(name).VEC.kass,res.(name).VEC.kass_SE,green);
    ylabel('nRV', 'Interpreter','latex')

    % plot nSI
    yyaxis right
    plotErrorbar(times,res.(name).nSI.kass,res.(name).nSI.kass_SE,viol);
    % plotErrorbar(times,res.(name).EVC.kass,res.(name).EVC.kass_SE,viol);
    ylabel('$n\Psi$', 'Interpreter','latex')

    xline(0,'r--')
    xlim([times(1),times(end)])
    xlabel('Time [ms]', 'Interpreter','latex')
    


    % plot nSI
    subplot(2,3,5)
    set(gca,'colororder',[viol;skyblue])
    yyaxis left
    plotErrorbar(times,res.(name).nSI.kass,res.(name).nSI.kass_SE,viol);
    % plotErrorbar(times,res.(name).EVC.kass,res.(name).EVC.kass_SE,viol);
    ylabel('$n\Psi$', 'Interpreter','latex')

    % plot CV2
    yyaxis right
    plotErrorbar(times,res.(name).cv2.cv2,res.(name).cv2.cv2_SE,skyblue);
    ylabel('Local CV$\textsuperscript{2}$', 'Interpreter','latex')
    ylim(areaParam(counterName).cv_lim)

    xline(0,'r--')
    xlim([times(1),times(end)])
    xlabel('Time [ms]', 'Interpreter','latex')


    sgtitle('Investigation of Auditory Regions of Birdsong for one Stim', 'Interpreter','latex')

end



