%按标签排序（单个表格单个工作表）
filename = 'E:\SEED\SEED_EEG\分开导入\实验融合\1.xlsx';
sheetName = 'Sheet1';
% 读取 Excel 数据
mergedData = readmatrix(filename, 'Sheet', sheetName);
% 获取第一列的值，之前需手动贴标签
firstColumn = mergedData(:, 1);
% 根据第一列的值将数据分成三部分（1、0、-1）
data_1 = mergedData(firstColumn == 1, :);
data_0 = mergedData(firstColumn == 0, :);
data_minus_1 = mergedData(firstColumn == -1, :);
% 合并这些部分
finalData = [data_1; data_0; data_minus_1];
% 将合并后的数据写入新 Excel 文件
% 指定 Excel 文件的路径和工作表名称
filename = 'E:\SEED\SEED_EEG\分开导入\实验融合\1.xlsx';
sheetName = 'Sheet1';
% 写入数据到 Excel
writematrix(finalData, filename, 'Sheet', sheetName);





% 所有表格所有工作表
sheetNames = {'Sheet1', 'Sheet2', 'Sheet3', 'Sheet4', 'Sheet5', 'Sheet6'};
% 定义要处理的Excel文件的编号
excelNumbers = 1:8;
for j = excelNumbers
    for i = 1:length(sheetNames)
        mergedData = readmatrix(['E:\SEED\SEED_EEG\分开导入\实验融合\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i});
        firstColumn = mergedData(:, 1);
        data_1 = mergedData(firstColumn == 1, :);
        data_0 = mergedData(firstColumn == 0, :);
        data_minus_1 = mergedData(firstColumn == -1, :);
        finalData = [data_1; data_0; data_minus_1];
        writematrix(finalData, ['E:\SEED\SEED_EEG\分开导入\实验融合\' num2str(j) '.xlsx'], 'Sheet', sheetNames{i});
    end
end

