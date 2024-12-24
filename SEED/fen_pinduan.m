% 设置输入和输出文件夹路径
inputFolder = 'E:\DEAP\sub\sub1\分好traila后有通道\';
outputFolder = 'E:\DEAP\sub\sub1\分频段\5';

% 确保输出文件夹存在，如果不存在则创建
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% 获取输入文件夹中的所有文件
fileList = dir(fullfile(inputFolder, '*.set'));
% 循环处理每个文件
for i = 1:numel(fileList)
    % 构建输入文件的完整路径
    inputFilePath = fullfile(inputFolder, fileList(i).name);
    % 加载 EEG 数据集
    EEG = pop_loadset('filename', fileList(i).name, 'filepath', inputFolder);
    % (delta: 1 -3Hz , theta: 4- 7Hz , alpha: 8- 13Hz , beta: 14 -30Hz，gamma : 31 -45Hz )
    EEG.data = eegfilt(EEG.data, EEG.srate,31,45); 
    % 构建输出文件的完整路径
    outputFilePath = fullfile(outputFolder, fileList(i).name);    
    % 保存处理后的 EEG 数据集到指定文件夹
    EEG = pop_saveset(EEG, 'filename', fileList(i).name, 'filepath', outputFolder);
end
