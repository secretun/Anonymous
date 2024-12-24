function [result,num_seg] = DataSeg(data,seg_len,overlap,ica,max_step,lrate,add_noise,p)
[row,col] = size(data);
% addnoise -- 0(default), 1
if add_noise
    R = wgn(row,col,p);
    data = data + R;
end

result = [];
index = 0;
num_seg = 0;
Fs = 250;
Fbp = [0.1,45];
N = 4;

% data = ft_preproc_bandpassfilter(data,Fs,Fbp,N);
% data = data - mean(data,2);
% data = data ./ std(data,0,2);
% ica -- 1(default),不对data做ICA处理
%        2,对全段data做ICA处理
%        3,对data_segment做ICA处理

switch ica
    case 1
        for i = 1:(seg_len-overlap):col
            if i+seg_len-1 <= col
                index = index + 1;
                num_seg = num_seg + 1;
                data_seg = data(:,i:i+seg_len-1);
                data_seg = data_seg - mean(data_seg,2);
                data_seg = data_seg ./ std(data_seg,0,2);
                data_filt = ft_preproc_bandpassfilter(data_seg,Fs,Fbp,N);
%                 data_filt = data_filt - mean(data_filt,2);
%                 data_filt = data_filt ./ std(data_filt,0,2);
%                 data_seg = data_seg - mean(data_seg,2);
%                 data_seg = data_seg ./ std(data_seg,0,2);
%                 result(:,:,index) = data_seg;
                result(:,:,index) = data_filt;
            end
        end
    case 2
        data = data - mean(data,2);
        data = data ./ std(data,0,2);
        [~,data_ica] = Infomax_ICA(data,max_step,lrate);
        for i = 1:(seg_len-overlap):col
            if i+seg_len-1 <= col
                index = index + 1;
                num_seg = num_seg + 1;
                data_seg = data_ica(:,i:i+seg_len-1);
                data_filt = ft_preproc_bandpassfilter(data_seg,Fs,Fbp,N);
%                 data_filt = data_filt - mean(data_filt,2);
%                 data_filt = data_filt ./ std(data_filt,0,2);
                result(:,:,index) = data_filt;
%                 data_seg = data_seg - mean(data_seg,2);
%                 data_seg = data_seg ./ std(data_seg,0,2);
%                 result(:,:,index) = data_seg;
            end
        end
    case 3
        for i = 1:(seg_len-overlap):col
            if i+seg_len-1 <= col
                index = index + 1;
                num_seg = num_seg + 1;
                data_seg = data(:,i:i+seg_len-1);
                data_seg = data_seg - mean(data_seg,2);
                data_seg = data_seg ./ std(data_seg,0,2);
                [~,data_ica] = Infomax_ICA(data_seg,max_step,lrate);
                data_filt = ft_preproc_bandpassfilter(data_ica,Fs,Fbp,N);
%                 data_filt = data_filt - mean(data_filt,2);
%                 data_filt = data_filt ./ std(data_filt,0,2);
%                 data_seg = data_seg - mean(data_seg,2);
%                 data_seg = data_seg ./ std(data_seg,0,2);
%                 result(:,:,index) = data_seg;
                result(:,:,index) = data_filt;
            end
        end
%     case 4
%         [data_iva,~] = ivabss(data,250,max_step,[],lrate,[]);
%         for i = 1:(seg_len-overlap):col
%             if i+seg_len-1 <= col
%                 index = index + 1;
%                 num_seg = num_seg + 1;
%                 data_seg = data_iva(:,i:i+seg_len-1);
%                 data_filt = ft_preproc_bandpassfilter(data_seg,Fs,Fbp,N);
%                 data_filt = data_filt - mean(data_filt,2);
%                 data_filt = data_filt ./ std(data_filt,0,2);
% %                 data_seg = data_seg - mean(data_seg,2);
% %                 data_seg = data_seg ./ std(data_seg,0,2);
% %                 result(:,:,index) = data_seg;
%                 result(:,:,index) = data_filt;
%             end
%         end
%     case 5
%         for i = 1:(seg_len-overlap):col
%             if i+seg_len-1 <= col
%                 index = index + 1;
%                 num_seg = num_seg + 1;
%                 data_seg = data(:,i:i+seg_len-1);
%                 [data_seg,~] = ivabss(data_seg,250,max_step,[],lrate,[]);
%                 data_filt = ft_preproc_bandpassfilter(data_seg,Fs,Fbp,N);
%                 data_filt = data_filt - mean(data_filt,2);
%                 data_filt = data_filt ./ std(data_filt,0,2);
% %                 data_seg = data_seg - mean(data_seg,2);
% %                 data_seg = data_seg ./ std(data_seg,0,2);
% %                 result(:,:,index) = data_seg;
%                 result(:,:,index) = data_filt;
%             end
%         end
end
end