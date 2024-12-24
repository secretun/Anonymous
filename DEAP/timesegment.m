clc; clear;
baseDir = 'E:\DEAP\sub\sub9\分频段\';
subDirs = {'1', '2', '3', '4', '5'};
if ~exist('E:\DEAP\sub\sub9\分时间段_频段\', 'dir')
    mkdir('E:\DEAP\sub\sub9\分时间段_频段\');
end
for s = 1:length(subDirs)
    if ~exist(fullfile('E:\DEAP\sub\sub9\分时间段_频段\', subDirs{s}), 'dir')
        mkdir(fullfile('E:\DEAP\sub\sub9\分时间段_频段\', subDirs{s}));
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
            save_path = fullfile('E:\DEAP\sub\sub9\分时间段_频段\', subDirs{s}, save_name);
            pop_saveset(EEG_segmented, 'filename', save_name, 'filepath', fullfile('E:\DEAP\sub\sub9\分时间段_频段\', subDirs{s}));
        end
    end
end