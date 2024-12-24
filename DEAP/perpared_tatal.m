clc;clear;
inputPath = 'E:\DEAP\data\维度变化后\s20.set';
outputFolder = 'E:\DEAP\sub\sub20\分好trail无通道';
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
clc;clear;
setDir = 'E:\DEAP\sub\sub20\分好trail无通道';
set1Dir = 'E:\DEAP\sub\sub20\分好trail有通道';
chanFilePath = 'E:\DEAP\chan\channel_32.mat';
if ~exist(set1Dir, 'dir')
    mkdir(set1Dir);
end
chanData = load(chanFilePath);
chanlocs = chanData.EEG.chanlocs;
setFiles = dir(fullfile(setDir, '*.set'));
eeglab;
for i = 1:length(setFiles)
    setFilePath = fullfile(setDir, setFiles(i).name);
    EEG = pop_loadset('filename', setFilePath);
    EEG.chanlocs = chanlocs;
    EEG = eeg_checkset(EEG);
    [~, name, ~] = fileparts(setFiles(i).name);
    newSetFilePath = fullfile(set1Dir, [name, '.set']);
    EEG = pop_saveset(EEG, 'filename', newSetFilePath);
end
clc;clear;
inputFolder = 'E:\DEAP\sub\sub20\分好trail有通道\';
freqBands = {[1 3], [4 7], [8 13], [14 30], [31 45]};
outputFolders = {'E:\DEAP\sub\sub20\分频段\1', 'E:\DEAP\sub\sub20\分频段\2', 'E:\DEAP\sub\sub20\分频段\3', 'E:\DEAP\sub\sub20\分频段\4', 'E:\DEAP\sub\sub20\分频段\5'};
fileList = dir(fullfile(inputFolder, '*.set'));
for i = 1:numel(fileList)
    inputFilePath = fullfile(inputFolder, fileList(i).name);
    EEG = pop_loadset('filename', fileList(i).name, 'filepath', inputFolder);
    for j = 1:length(freqBands)
        if ~exist(outputFolders{j}, 'dir')
            mkdir(outputFolders{j});
        end
        EEG_filtered = EEG;
        EEG_filtered.data = eegfilt(EEG.data, EEG.srate, freqBands{j}(1), freqBands{j}(2));
        outputFilePath = fullfile(outputFolders{j}, fileList(i).name);
        EEG_filtered = pop_saveset(EEG_filtered, 'filename', fileList(i).name, 'filepath', outputFolders{j});
    end
end
clc; clear;
baseDir = 'E:\DEAP\sub\sub20\分频段\';
subDirs = {'1', '2', '3', '4', '5'};
if ~exist('E:\DEAP\sub\sub20\分时间段_频段\', 'dir')
    mkdir('E:\DEAP\sub\sub20\分时间段_频段\');
end
for s = 1:length(subDirs)
    if ~exist(fullfile('E:\DEAP\sub\sub20\分时间段_频段\', subDirs{s}), 'dir')
        mkdir(fullfile('E:\DEAP\sub\sub20\分时间段_频段\', subDirs{s}));
    end
end
for s = 1:length(subDirs)
    set_folder = fullfile(baseDir, subDirs{s});
    set_files = dir(fullfile(set_folder, '*.set'));
    for i = 1:length(set_files)
        set_path = fullfile(set_folder, set_files(i).name);
        EEG = pop_loadset('filename', set_files(i).name, 'filepath', set_folder);
        segment_duration = 20;
        data_points = EEG.pnts;
        sampling_rate = EEG.srate;
        data_points_per_segment = segment_duration * sampling_rate;
        num_segments = floor(data_points / data_points_per_segment);
        for j = 1:num_segments
            segment_start = (j - 1) * data_points_per_segment + 1;
            segment_end = j * data_points_per_segment;
            segmented_data = EEG.data(:, segment_start:segment_end);
            EEG_segmented = EEG;
            EEG_segmented.data = segmented_data;
            EEG_segmented.pnts = size(segmented_data, 2);
            [~, set_name, set_ext] = fileparts(set_files(i).name);
            save_name = [set_name '_segment_' num2str(j) set_ext];
            save_path = fullfile('E:\DEAP\sub\sub20\分时间段_频段\', subDirs{s}, save_name);
            pop_saveset(EEG_segmented, 'filename', save_name, 'filepath', fullfile('E:\DEAP\sub\sub20\分时间段_频段\', subDirs{s}));
        end
    end
end