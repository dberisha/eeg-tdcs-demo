% prepData

% This script will help convert .xdf into .mat file by doing the following

% 1. Loading XDF
% 2. Figure out eeg adn trigger
% 3. Truncate base on trigger into 7 data segments
% 4. Save out those 7 individual data setments into 7 mat files

% - DTN: 3/22/2019 (Destiny)

%% Clearing all variables and setting directories

ccc;

inDir = '/Users/destinyberisha/Desktop/eeg-tdcs-demo/data/raw/xdf';
outDir = '/Users/destinyberisha/Desktop/eeg-tdcs-demo/data/processed/epochedMat';

cd(inDir)

%% Look Up and Load XDF Files

% subID = input('What is the subject''s id?', s');

subData = dir('*.xdf');

[streams, fileheader] = load_xdf(subData.name);


%% Finding the correct data streams

streams{1}.info; % Getting information of the first cell, curly brackets special for cells

streams{1}; % Evaluate to see what this is 

% Markers are the triggers, series is the data of the EEG
% Always check which is which 


for i =1:length(streams)
    if strcmp(streams{i}.info.name,'BrainAmpSeries')
        eegInd=i;
    else
        trigInd=i;
    end
end


%% Extracting out the eeg data and info

eeg = (streams{eegInd}.time_series).';
eeg_time = streams{eegInd}.time_stamps;
nEEGsamples = numel(eeg_time);
fs_eeg = 1/mean(diff(streams{eegInd}.time_stamps));


%% Extracting the markers, trigger data info


markers = streams{trigInd}.time_series;
markers_time = streams{trigInd}.time_stamps;
nMarkers = numel(markers_time(:));
trigger = zeros(nMarkers,1);
for m=1:nMarkers
    trigger(m) = str2num(markers{m});
end

%% Mark off data segments for each trial, preparing for truncation



fprintf('\nMark the start of each trial block.\n.')
figure;
plot(markers_time, trigger,'.');
[starts,~] = ginput(7);

fprintf('\nMark the end of each trial block.\n')
figure;
plot(markers_time, trigger, '.');
[ends,~] = ginput(7);

close all;










