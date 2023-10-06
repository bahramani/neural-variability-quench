% This code will transform each data so that var_analysis.m can accept them
close all
clear
clc

refData = load('SNR.mat');

%% IpN
name = 'IpN';
load("IpN.mat")

tmp = data.spikeTimes';
rasterAll = cell(size(tmp));

for i = 1:length(tmp)
    rasterAll{i} = cell(length(tmp{i}(1,:)),1);
    for j = 1:length(rasterAll{i})
        rasterAll{i}{j} = tmp{i}(:,j);
        rasterAll{i}{j} = rasterAll{i}{j}(~isnan(rasterAll{i}{j}));
        rasterAll{i}{j} = rasterAll{i}{j} .* 1000;
    end
end

raster = rasterAll(1:270);
save(['Data_Trans/',name,'1.mat'], 'raster')

raster = rasterAll(271:372);
save(['Data_Trans/',name,'2.mat'], 'raster')

%% Allen VISam

load("VISam.mat")
name = 'VIsam';

[~, numTrial, numNeuron] = size(visam);
raster = cell(numNeuron,1);
for i = 1:numNeuron
    raster{i} = cell(numTrial,1);
    for j = 1:numTrial
        raster{i}{j} = find(visam(:,j,i))-500;
    end
end
save(['Data_Trans/',name,'.mat'], 'raster')

%% Allen VISp

load("VISp.mat")
name = 'VISp';

[~, numTrial, numNeuron] = size(visp);
raster = cell(numNeuron,1);
for i = 1:numNeuron
    raster{i} = cell(numTrial,1);
    for j = 1:numTrial
        raster{i}{j} = find(visp(:,j,i))-500;
    end
end
save(['Data_Trans/',name,'.mat'], 'raster')

%% Allen MGd

load("MGd.mat")
name = 'MGd';

[~, numTrial, numNeuron] = size(MGd);
raster = cell(numNeuron,1);
for i = 1:numNeuron
    raster{i} = cell(numTrial,1);
    for j = 1:numTrial
        raster{i}{j} = find(MGd(:,j,i))-500;
    end
end
save(['Data_Trans/',name,'.mat'], 'raster')

%%








