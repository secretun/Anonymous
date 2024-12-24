%按照标签进行重新排序
filename = 'E:\SEED\microstate\分开导入\实验融合\fea1.xlsx';
%sheetNames = {'Sheet1', 'Sheet2', 'Sheet3', 'Sheet4', 'Sheet5', 'Sheet6'};
sheetNames = {'Sheet6'};
% Loop through each sheet
for i = 1:length(sheetNames)
    % Read data from the current sheet
    mergedData = readmatrix(filename, 'Sheet', sheetNames{i});
    % Sort data
    firstColumn = mergedData(:, 1);
    data_1 = mergedData(firstColumn == 1, :);
    data_0 = mergedData(firstColumn == 0, :);
    data_minus_1 = mergedData(firstColumn == -1, :);
    finalData = [data_1; data_0; data_minus_1];
    % Write sorted data back to the same sheet
    writematrix(finalData, filename, 'Sheet', sheetNames{i});
end




%特征融合 需手动删除Sheet6中4个标签列
filePath = 'E:\SEED\microstate\分开导入\第三次实验\15.xlsx';
excelFile = xlsread(filePath);
mergedData = [];
for sheetIndex = 1:5
    currentData = xlsread(filePath, sheetIndex);
    mergedData = [mergedData, currentData];
end 
xlswrite(filePath, mergedData, 'Sheet6');