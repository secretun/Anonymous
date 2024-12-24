function de = DEOf5Bands(data,step,overlap)
    % data : channel*length*segment_num
    delta = [1,3];
    theta = [4,7];
    alpha = [8,13];
    beta = [14,30];
    gamma = [31,50];
    Fs = 200; N = 4;
    [~,len,seg_num] = size(data);
    for i = 1:seg_num
        data_delta(:,:,i) = ft_preproc_bandpassfilter(data(:,:,i),Fs,delta,N); 
        data_theta(:,:,i) = ft_preproc_bandpassfilter(data(:,:,i),Fs,theta,N); 
        data_alpha(:,:,i) = ft_preproc_bandpassfilter(data(:,:,i),Fs,alpha,N); 
        data_beta(:,:,i)  = ft_preproc_bandpassfilter(data(:,:,i),Fs,beta,N); 
        data_gamma(:,:,i) = ft_preproc_bandpassfilter(data(:,:,i),Fs,gamma,N); 
    end
    count = 0;
    for i = 1:seg_num
        for j = 1:(step-overlap):len
            if j+step-1 <= len
                count = count + 1;  
                seg_delta = data_delta(:,j:(j+step-1),i);
                seg_theta = data_theta(:,j:(j+step-1),i);
                seg_alpha = data_alpha(:,j:(j+step-1),i);
                seg_beta = data_beta(:,j:(j+step-1),i);
                seg_gamma = data_gamma(:,j:(j+step-1),i);
                
                std_delta = std(seg_delta,0,2);
                std_theta = std(seg_theta,0,2);
                std_alpha = std(seg_alpha,0,2);
                std_beta  = std(seg_beta,0,2);
                std_gamma = std(seg_gamma,0,2);
                
                de_delta(:,count) = 1/2 * (log(2.*pi.*std_delta.^2) + 1);
                de_theta(:,count) = 1/2 * (log(2.*pi.*std_theta.^2) + 1);
                de_alpha(:,count) = 1/2 * (log(2.*pi.*std_alpha.^2) + 1);
                de_beta(:,count) = 1/2 * (log(2.*pi.*std_beta.^2) + 1);
                de_gamma(:,count) = 1/2 * (log(2.*pi.*std_gamma.^2) + 1);
            end
        end
    end
    de = [de_delta;de_theta;de_alpha;de_beta;de_gamma];

end