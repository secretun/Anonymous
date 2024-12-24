
inputFolder = 'E:\DEAP\sub\sub9\分好trail有通道\';
freqBands = {[1 3], [4 7], [8 13], [14 30], [31 45]};
outputFolders = {'E:\DEAP\sub\sub9\分频段\1', 'E:\DEAP\sub\sub9\分频段\2', 'E:\DEAP\sub\sub9\分频段\3', 'E:\DEAP\sub\sub9\分频段\4', 'E:\DEAP\sub\sub9\分频段\5'};
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