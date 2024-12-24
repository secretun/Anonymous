% 设置输入和输出文件夹路径
inputFolder = 'E:\SEED\SEED_EEG\官网SEED\Preprocessed_EEG\'; % 包含.mat文件的文件夹
outputFolder = 'E:\SEED\SEED_EEG\nochanol\15_20131105\'; % 用于保存.set文件的文件夹

% 指定要处理的.mat文件
fileToProcess = '15_20131105'; % 用实际文件名替换 'your_specific_file.mat'

% 构建.mat文件的完整路径
matFilePath = fullfile(inputFolder, fileToProcess);

% 加载.mat文件
loadeddata=load(matFilePath);


% 创建一个新的 EEG 结构
EEG = eeg_emptyset;

% 构建输出 .set 文件的完整路径
[~, baseFileName, ~] = fileparts(fileToProcess);

% 循环赋值给 EEG.data
for j = 1:15  % 假设有3个不同编号的 EEG 数据文件，你可以根据实际情况调整循环次数
    % 根据文件名的编号构建变量名，例如 eeg1、eeg2、eeg3
    eegVarName = ['zjy_eeg' num2str(j)];

    % 从工作区中获取数据并分配给 EEG.data
    EEG.data = eval(eegVarName); 

    % 构建输出 .set 文件的完整路径，例如 Data1.set、Data2.set、Data3.set
    setFilePath = fullfile(outputFolder, [baseFileName 'Data' num2str(j) '.set']);

    % 使用 EEGLAB 函数保存为 .set 文件
    pop_saveset(EEG, setFilePath);
    
    % 清除 EEG.data，以便下一次循环赋值
    EEG.data = [];
end



% 清除工作区中的变量
clear EEG;
