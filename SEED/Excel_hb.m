%将所有被试表格合并成一个表格（45x15,40）内存不足
% 设置表格文件夹路径
folderPath = 'E:\SEED\SEED_EEG\分频段表格\';

% 创建一个新的表格来保存合并的结果
combinedTable = [];

% 循环遍历表格1到15
for i = 1:15
    % 生成表格文件名，假设文件名是数字
    fileName = [num2str(i) '.xlsx'];
    
    % 读取当前表格
    currentTable = readtable(fullfile(folderPath, fileName));
    
    % 将当前表格添加到合并表格
    if i == 1
        combinedTable = currentTable;
    else
        combinedTable = [combinedTable; currentTable];
    end
end

% 将合并的表格保存为表格16
writetable(combinedTable, fullfile(folderPath, '表格16.xlsx'));
