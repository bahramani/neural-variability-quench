%% This code generates data for Plot_results.m (in case data is not already be available in /Analyzed folder)
% Citation: Saleh F, Fakharian M & Ghazizadeh A "Stimulus presentation can enhance spiking irregularity across subcortical and cortical regions." PLoS Comp Biol, 2022

addpath(genpath('C:\Users\user\Desktop\Amirreza\Code_and_Data\lib'))
%%
clc;
clear;

names = {'stim_all', 'stim_10', 'stim_11', 'stim_12', 'stim_16', 'stim_44', 'stim_45', 'stim_46', 'stim_50', 'stim_51', 'stim_52'};

for counter_name = 2:numel(names)
    
    % params and region selection
    
    clearvars -except counter_name names
    name = names{counter_name};
    disp(name)
    
    load(['Datasets/',name,'.mat'])
    
    switch(name)

        case {'MGd', 'VIsam', 'VISp'}
            % pre and post time
            T = [-200,300];

            align = 0;

            params.kass_params.K = 20;

            % params
            params.t_start = -500;
            params.t_end = 750;
            params.bin_width = 500;
            params.slide = 50;


        case{'stim_all', 'stim_10', 'stim_11', 'stim_12', 'stim_16', 'stim_44', 'stim_45', 'stim_46', 'stim_50', 'stim_51', 'stim_52'}
            % pre and post time
            T = [-200,500];

            align = 0;

            params.kass_params.K = 20;

            % params
            params.t_start = -200;
            params.t_end = 1000;
            params.bin_width = 300;
            params.slide = 50;

        case('IpN1')
            % pre and post time
            T = [-200,250];

            align = 0;

            params.kass_params.K = 20;

            % params
            params.t_start = -500;
            params.t_end = 1500;
            params.bin_width = 300;
            params.slide = 50;

        case('IpN2')
            % pre and post time
            T = [-100,350];

            align = 0;

            params.kass_params.K = 20;

            % params
            params.t_start = -200;
            params.t_end = 800;
            params.bin_width = 300;
            params.slide = 50;
        
        case('vlPFC')
            % pre and post time
            T = [-100,50];
            
            align = 0;
            
            params.kass_params.K = 20;
            
            % params
            params.t_start = -200;
            params.t_end = 800;
            params.bin_width = 200;
            params.slide = 50;
            
        case('SNR')
            % pre and post time
            T = [-100,50];
            
            params.kass_params.K = 20;
            
            align = 0;
            
            % params
            params.t_start = -200;
            params.t_end = 800;
            params.bin_width = 200;
            params.slide = 50;
            
    end
    % Analysis
    raster = raster(:);
    n = numel(raster);
    t_end = params.t_end;
    t_start = params.t_start;
    bin_width = params.bin_width;
    slide = params.slide;
    min_spike = 3;
    
    % Removing low rate neurons
    % find average spike count of each neuron
%     valid_index = ones(n,1);
%     
%     for i = 1:n
%         raster_local = raster{i};
%         spike_mat = cell2mat(cellfun(@(x) ...
%             histcounts(x,t_start:bin_width:t_end),...
%             raster_local,'UniformOutput',false));
%         if mean(spike_mat(:)) < min_spike
%             valid_index(i) = 0;
%         end
%     end
%     
%     raster = raster(valid_index == 1);
%     n = length(raster);
%     disp(sum(valid_index)/length(valid_index)*1e2);
    
    % rate
    data = cell(n,1);
    
    for i = 1:n
        raster_local = raster{i};
        spike_mat = cell2mat(cellfun(@(x) histcounts(x,t_start:t_end),...
            raster_local,'UniformOutput',false));
        data{i} = spike_mat;
    end
    
    bin_size = 10;
    
    fr = [];
    
    for i = 1:n
        perccount(i,numel(raster))
        fr_temp = mean(data{i});
        fr(i,:) = movmean(fr_temp,bin_size,'Endpoints','discard')*1e3;
    end
    
    rate = mean(downsample(fr,2));
    t_rate = t_start:t_end-bin_size;
    
    % FFvT
    bin_ratio = 0:0.01:1;
    
    parfor i = 1:numel(raster)
        perccount(i,numel(raster))
        
        raster_local = raster{i};
        
        [~,ypre(i,:)] = plot_slope_resampled(raster_local,T(1),bin_width,bin_ratio);
        [X,ypost(i,:)] = plot_slope_resampled(raster_local,T(2),bin_width,bin_ratio);
    end
    
    ind = any(isnan(ypre(:,2:end)),2) | any(isnan(ypost(:,2:end)),2);
    
    ypre(ind,:) = [];
    ypost(ind,:) = [];
    
    FFpre_mean = mean(ypre);
    FFpost_mean = mean(ypost);
    FFpre_SE = nanstd(ypre)/sqrt(length(ypre));
    FFpost_SE = nanstd(ypost)/sqrt(length(ypost));
    
    slope_pre = zeros(length(ypre),1);
    intercept_pre = zeros(length(ypre),1);
    slope_post = zeros(length(ypre),1);
    intercept_post = zeros(length(ypre),1);
    
    tbin_part = bin_ratio(71:end);
    ypre_part = ypre(:,71:end);
    ypost_part = ypost(:,71:end);
    
    parfor i = 1:size(ypre,1)
        fit_pre = fitlm(tbin_part,ypre_part(i,:));
        fit_post = fitlm(tbin_part,ypost_part(i,:));
        slope_pre(i) = fit_pre.Coefficients.Estimate(2);
        slope_post(i) = fit_post.Coefficients.Estimate(2);
        intercept_pre(i) = fit_pre.Coefficients.Estimate(1);
        intercept_post(i) = fit_post.Coefficients.Estimate(1);
    end
    
    inter_pre = mean(intercept_pre);
    inter_pre_SE = std(intercept_pre)/sqrt(length(ypre));
    sl_pre = mean(slope_pre);
    sl_pre_SE = std(slope_pre)/sqrt(length(ypre));
    
    inter_post = mean(intercept_post);
    inter_post_SE = std(intercept_post)/sqrt(length(ypre));
    sl_post = mean(slope_post);
    sl_post_SE = std(slope_post)/sqrt(length(ypre));

    cv2 = var_cv(raster,params);
    params.kass_params.cv2_all_neu = cv2.cv2_all_neu;
    results = var_decom(raster,params);
    results.cv2 = cv2;
    results.FFvT.pre.FF_mean = FFpre_mean;
    results.FFvT.pre.FF_SE = FFpre_SE;
    results.FFvT.post.FF_mean = FFpost_mean;
    results.FFvT.post.FF_SE = FFpost_SE;
    results.FFvT.T = bin_ratio*bin_width;
    results.FFvT.sl_pre = sl_pre;
    results.FFvT.sl_pre_SE = sl_pre_SE;
    results.FFvT.sl_post = sl_post;
    results.FFvT.sl_post_SE = sl_post_SE;
    results.FFvT.inter_pre = inter_pre;
    results.FFvT.inter_pre_SE = inter_pre_SE;
    results.FFvT.inter_post = inter_post;
    results.FFvT.inter_post_SE = inter_post_SE;
    results.rate.fr = rate;
    results.rate.t = t_rate;
    save(['Analyzed2/',name,'_Analyzed.mat'],'results')
end