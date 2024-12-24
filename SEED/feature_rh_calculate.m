%导入每个被试的5个频段数据（3x15x5)个数据，与feature_fusion1连用
tic
clear;  clc;
% start EEGLAB to load all dependent paths
eeglab
% 指定基本路径和文件名
baseDir = 'E:\SEED\SEED_EEG\分频段数据\1_20131030\';
filenames = dir(baseDir);
%all_feature_data = zeros(5,15,40);
% 循环遍历文件名
for i = 3:length(filenames)
    % 构建当前文件夹的完整路径
    filename = dir([baseDir,filenames(i).name,filesep,'*.set']);
        % 循环遍历数据集
    for j= 1:15
        % 构建完整文件路径
        trial_name = filename(j).name;
        trial_root = filename(j).folder;
        % 示例：加载数据
        EEG = pop_loadset('filename',trial_name, 'filepath', trial_root);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
    end
end
  %% 3.3.2 Select data for microstate analysis  计算GFP
    [EEG, ALLEEG] = pop_micro_selectdata( EEG, ALLEEG, 'datatype', 'spontaneous',...
    'avgref', 1, ...
    'normalise', 1, ...
    'MinPeakDist', 10, ...
    'Npeaks', 1000, ...
    'GFPthresh', 1, ...
    'dataset_idx', 1:225 );
    % store data in a new EEG structure
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw % updates EEGLAB datasets
    %% 3.4 Microstate segmentation  对GFP进行聚类
    % select the "GFPpeak" dataset and make it the active set
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 225,'retrieve',226,'study',0);
    eeglab redraw
    % Perform the microstate segmentation 执行微状态分割 
    EEG = pop_micro_segment( EEG, 'algorithm', 'modkmeans', ...
    'sorting', 'Global explained variance', ...
    'Nmicrostates', 2:8, ...
    'verbose', 1, ...
    'normalise', 0, ...
    'Nrepetitions', 50, ...
    'max_iterations', 1000, ...
    'threshold', 1e-06, ...
    'fitmeas', 'CV',...
    'optimised',1);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    %% 3.5 Review and select microstate segmentation
    %% 3.5.1 Plot microstate prototype topographies
    figure;MicroPlotTopo( EEG, 'plot_range', [] );
    %% 3.5.2 Select active number of microstates
    EEG = pop_micro_selectNmicro( EEG, 'Nmicro', 5 );
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
     occurence = cell(225,1);
    duration = cell(225,1);
    coverage = cell(225,1);
    MGfp = cell(225,1);
    data = cell(225,1);
    TP_all = cell(225,1);
    fit = cell(225,1);
    prototypes = EEG.microstate.prototypes;
    for k = 1:225
        fprintf('Importing prototypes and backfitting for dataset %i\n',k)
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',k,'study',0);
    %% Import microstate prototypes from other dataset to the datasets that should be back-fitted
    % 将其他数据集的微状态原型导入到需要进行回配的数据集  
    % note that dataset number 5 is the GFPpeaks dataset with the microstate
    % prototypes
        EEG = pop_micro_import_proto( EEG, ALLEEG,226);
    %% 3.6 Back-fit microstates on EEG
        EEG = pop_micro_fit( EEG, 'polarity', 0 );
     
    %% 3.7 Temporally smooth microstates labels
        EEG = pop_micro_smooth( EEG, 'label_type', 'backfit', ...
        'smooth_type', 'reject segments', ...
        'minTime', 30, ...
        'polarity', 0 );
    %% 3.9 Calculate microstate statistics
        EEG = pop_micro_stats( EEG, 'label_type', 'backfit', ...
        'polarity', 0 );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        occurence{k,1} = EEG.microstate.stats.Occurence;
        duration{k,1} = EEG.microstate.stats.Duration;
        coverage{k,1} = EEG.microstate.stats.Coverage;
        fitlabels = EEG.microstate.fit.labels;
        MGfp{k,1} =EEG.microstate.stats.Gfp;
        TP = EEG.microstate.stats.TP;
        GEVall(k,:) = EEG.microstate.stats.GEV;
        TP_all{k,1} = TP; % 按列排序
    %     data{i,1} = [occurence, duration, coverage];
        fit{k,1} = fitlabels;
    end
    toc
     disp(['程序运行时间', num2str(toc)]) ;
    %% 3.8 Illustrating microstate segmentation
    %Plotting GFP of active microstates for the first 1500 ms for subject 1.
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',1,'study',0);
    figure;MicroPlotSegments( EEG, 'label_type', 'backfit', ...
    'plotsegnos', 'first', 'plot_time', [1 100000], 'plottopos', 1 );
    eeglab redraw
     %合并数据
    features_data = feature_fusion1(coverage, duration, occurence,TP_all);
    
   %m=i-2;
   % all_feature_data(m,:,:) = features_data;
     xlswrite('E:\SEED\SEED_EEG\分频段表格\16.xlsx', features_data,1, 'B1'); 