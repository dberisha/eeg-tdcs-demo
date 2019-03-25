% prepData

% This script will help convert .xdf into .mat file by doing the following

% 1. Loading XDF
% 2. Figure out eeg adn trigger
% 3. Truncate base on trigger into 7 data segments
% 4. Save out those 7 individual data setments into 7 mat files

% - DTN: 3/22/2019 (Destiny)


% Notes by Destiny:

% First we are loading in the streams, and extracting the MARKERS from the
% DAQ. Then we are manually saving (and checking) the blocks that are
% aligned with these MARKERS. Then we are extracting the tiemstamps from
% the EEG data, and matching those up. Then matching up the blocks of the
% EEG timestamps to the actual EEG data. 

% Remember the timestamps are coming from the universal clock. 

%% Clearing all variables and setting directories


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
%% Calculate the epoch - segmenting trigger info for trial blocks

segEndInd=zeros(1,7);
segEngTime=zeros(1,7);
epochedTriggerCell = {};

for k=1:length(starts)
    b=1;
    for i=1:length(markers_time)
        if markers_time(i)>=starts(k) && markers_time(i)<ends(k)
        epochedTriggerCell{k}(b)= markers_time(i);
        segEngTime(k)=markers_time(i);
        segEndInd(k)=i;
        b=b+1;
        end
    end
end


 %% Sanity Check
 
 for i=1:length(epochedTriggerCell)
     A(i)=size(epochedTriggerCell{i},2);
 end
 
 prompt=sprintf('\nAre the num of the sizes of the segmented chunks (%d) equal to the original size (%d)? (y/n). \n', sum(A), size(markers_time,2));
 check = input(prompt, 's');
 
 if strcmp(check,'y')
     fprintf('\nYou''re good to go.')
 else
     fprintf('\nYou probably marked off data segments incorrectly in the figure, go back and try again.\n')
 end
 
%% Epoching EEG based on segmented trial blocks


epochedEEGCell={};
epochedEEGTimeCell={};
sFlashTime=1; % adding 1 second to the final cut of the trial block to be safe

for i=1:length(epochedTriggerCell)
    sTime=epochedTriggerCell{i}(1);
    eTime=epochedTriggerCell{i}(end) + sFlashTime;
    [~,sInd]=min(abs(eeg_time-sTime));
    [~,eInd]=min(abs(eeg_time-eTime));
    
    epochedEEGCell{i}=eeg(sInd:eInd,:);
    epochedEEGTimeCell{i}=eeg_time(sInd:eInd);
end

%% Splitting up the data and saving it out

subID='S01';
condStr='Anodal';

fprintf('\n Saving...Please be patient. \n');

cd(outDir)

for i=1:7
    
    epochedEEG = epochedEEGCell{i};
    epochedTrigger = epochedTriggerCell{i};
    epochedEEGTime =epochedEEGTimeCell{i};
    
    if i<=2
        endTag = '_preStim';
        bTag = sprintf('_b%02d.mat',i);
    elseif i>=3 && i<=5
        endTag = '_durStim';
        bTag = sprintf('_b%02d.mat',i);
    elseif i>=6
        endTag = '_posStim';
        bTag = sprintf('_b%02d.mat',i);
    end
    
    epochedFilename=[subID condStr endTag bTag];
    save(epochedFilename, 'epochedEEG', 'epochedTrigger','epochedEEGTime','fs_eeg', 'subID', 'condStr');
    clear epochedEEG epochedTrigger epochedEEGTime
end

% -7.3v
% add this tag if you need to save out more than 2GB

fprintf('\n Save Complete.\n')






