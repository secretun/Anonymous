%特征融合 需手动删除Sheet6中4个标签列
% 设置文件路径
filePath = 'E:\SEED\SEED_EEG\分开导入\第三次实验\8.xlsx';
% 打开Excel文件
excelFile = xlsread(filePath);
% 创建一个新的工作表，初始化为一个空矩阵
mergedData = [];
% 循环读取每个工作表
for sheetIndex = 1:5
    % 读取当前工作表的数据
    currentData = xlsread(filePath, sheetIndex);
    % 合并数据到新工作表，将数据逐行添加
    mergedData = [mergedData, currentData];
end
% 将合并后的数据写入一个新工作表
xlswrite(filePath, mergedData, 'Sheet6');









%合并被试3次实验特征值（竖）
sheetNames = {'Sheet1', 'Sheet2', 'Sheet3', 'Sheet4', 'Sheet5', 'Sheet6'};
% 定义要处理的Excel文件的编号
excelNumbers = 1:8;
for j = excelNumbers
    for i = 1:length(sheetNames)
 data1 = readtable(['E:\SEED\SEED_EEG\分开导入\第一次实验\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i});
 data2 = readtable(['E:\SEED\SEED_EEG\分开导入\第二次实验\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i});
 data3 = readtable(['E:\SEED\SEED_EEG\分开导入\第三次实验\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i});
        mergedData = vertcat(data1, data2, data3);
        writetable(mergedData, ['E:\SEED\SEED_EEG\分开导入\实验融合\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i}, 'WriteVariableNames', false);
    end
end
