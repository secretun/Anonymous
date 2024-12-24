%分trail
clc;clear;
inputPath = 'E:\DEAP\data\维度变化后\s11.set';
outputFolder = 'E:\DEAP\sub\sub dls\分好trail无通道';
EEG = pop_loadset(inputPath);
[numChannels, numTrials,numSamples] = size(EEG.data);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i = 1:numTrials
    newEEG = EEG;
    newEEG.data = EEG.data(:, i, :);
    newEEG.data = squeeze(newEEG.data);
    [~, fileName, ~] = fileparts(inputPath);
    newFileName = sprintf('%s_trial%d.set', fileName, i);
    outputPath = fullfile(outputFolder, newFileName);
    pop_saveset(newEEG, 'filename', outputPath);
end