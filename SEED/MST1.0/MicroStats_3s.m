%MICROSTATS Calculates microstate statistics.
%
% Usage:
%   >> Mstats = MicroStats(X, A, L)
%   >> Mstats = MicroStats(X, A, L,polarity,fs)
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
%  Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
%  Inputs:
%   X - EEG (channels x samples (x trials)).
%   A - Spatial distribution of microstate prototypes (channels x  K).微观状态空间聚类
%   L - Label of the most active microstate at each timepoint (trials x
%       time).每个时间点上最活跃的微状态的标签
%
%  Optional input:
%   polarity - Account for polarity when fitting Typically off for
%              spontaneous EEG and on for ERP data (default = 0).
%   fs       - Sampling frequency of EEG (default = 1).
%
%  Outputs:
%  Mstats - Structure of microstate parameters per trial:每次试验的微观状态参数结构
%   .Gfp        - Global field power
%   .Occurence  - Occurence of a microstate per s 每秒出现的微状态
%   .Duration   - Average duration of a microstate 微状态的平均持续时间
%   .Coverage   - % of time occupied by a microstate 被微观状态占用的时间
%   .GEV        - Global Explained Variance of microstate
%   .MspatCorr  - Spatial Correlation between template maps and microstates
%   .TP         - transition probabilities
%   .seq        - Sequence of reoccurrence of MS ((trials x)  time).
%   .msFirstTF  - First occurence of a microstate (similar to use like seq)
%                 ((trials x)  time)
%   .polarity   - see inputs.
%
%  Mstats.avgs - Structure of microstate parameters mean / and std over
%                trials.试验中微观状态参数均值和标准差的结构。
%
% Authors:
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Z黵ich, Psychologisches Institut, Methoden der
% Plastizit鋞sforschung.
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% September 2017.

% Copyright (C) 2017 Andreas Pedroni, andreas.pedroni@uzh.ch.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function Mstats = MicroStats_3s(X,A,L,polarity,fs)
%% Error check and initialisation 错误检查和初始化
if nargin < 5
    fs = 1;
elseif nargin < 3
    polarity = 0;
elseif nargin < 3
    help MicroStats;
    return;
end

[C,N,Ntrials] = size(X); % 也就说它把二维矩阵当作第三维为1的三维矩阵，这也如同我们把n维列向量当作n×1的矩阵一样,即Ntrials为1
K = size(A,2);  % 返回的是A矩阵的列数，即微状态原型的数目


%% Spatial correlation between microstate prototypes and EEG 微状态原型与脑电图的空间相关性
% force average reference
X = X - repmat(mean(X,1),[C,1,1]);  % mean(X,1)表示返回返回X矩阵每列的平均值1*samples，repmat(mean(X,1),[C,1,1])矩阵变为62*samples

%Check if the data has more than one trial or not and reshape if necessary 检查数据是否有不止一次的trail，必要时重新整形
if Ntrials > 1 % for epoched data
    X = squeeze(reshape(X, C, N*Ntrials));  % squeeze函数用于删除矩阵中的单一维，随机产生一个1x2x3的矩阵A，然后squeeze（A）将返回一个2x3的矩阵，将第一维却掉(因为第一位大小为1)
end

% Normalise EEG and maps (average reference and gfp = 1 for EEG)
Xnrm = X ./ repmat(std(X,1), C, 1); % already have average reference 已经有了平均参考 std(X,1)计算X的标准差，1表示除以N
A_nrm = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);

% Global map dissimilarity
GMD = nan(K,N);
for k = 1:K
    GMD(k,:) = sqrt(mean( (Xnrm - repmat(A_nrm(:,k),1,N)).^2 ));
end

% Account for polarity (recommended 0 for spontaneous EEG)
if polarity == 0
    GMDinvpol = nan(K,N);
    for k = 1:K
        GMDinvpol(k,:) = sqrt(mean( (Xnrm - repmat(-A_nrm(:,k),1,size(Xnrm,2))).^2));
    end
    idx = GMDinvpol < GMD;
    GMD(idx) = GMDinvpol(idx);
end

% Calculate the spatial correlation between microstate prototypes and EEG 计算微状态原型和脑电图之间的空间相关性
SpatCorr = 1 - (GMD.^2)./2;

% APedroni added this 13.12.2017
GEVtotal = sum((SpatCorr(sub2ind(size(SpatCorr),L,1:size(L,2))).*squeeze(std(X)).^2) ./sum(squeeze(std(X)).^2));

SpatCorr = squeeze(reshape(SpatCorr,K,N,Ntrials));

%% Sequentialize
% take into account the sequence of occurence of microstates. This makes
% only sense in ERP data (if it makes sense) 考虑微观状态的出现顺序。这只在ERP数据中有意义(如果有意义的话)
seq = nan(Ntrials, N);  % Ntrials为1，N为samples
msFirstTF = nan(Ntrials, N);
for trial = Ntrials
    first = 1;
    s = ones(K,1);  % 生成5*1的全1矩阵
    for n = 1:(N-1)
        if L(trial,n) == L(trial,n+1)  % L表示label
            seq(trial,n) = s(L(trial,n));
            msFirstTF(trial,n) = first;
        else
            seq(trial,n) = s(L(trial,n))  ;
            s(L(trial,n),1) = s(L(trial,n),1) + 1;
            first = n;
            msFirstTF(trial,n) = first;
        end
    end
    seq(trial,n+1) = seq(trial,n);
    msFirstTF(trial,n+1) = msFirstTF(trial,n);
end

%% 计算整个trail的MspatCorr和GEV
MspatCorr = nan(Ntrials,K);
GEV = nan(Ntrials,K);
GFP = squeeze(std(X));  % std(X)函数求解的是最常见的标准差,此时除以的是N-1，默认为按列，squeeze去除size为1的维度
if Ntrials>1
    GFP = GFP';
end
for k = 1:K
    for trial = 1:Ntrials
         % Average spatial correlation per Microstate 每个微状态的平均空间相关性
        MspatCorrTMP = SpatCorr(:,:,trial);
        MspatCorr(trial,k) = nanmean(MspatCorrTMP(k,L(trial,:)==k));
        
        % global explained variance Changed by Pedroni 3.1.2018
        GEV(trial,k) = sum( (GFP(trial,L(trial,:)==k) .* MspatCorrTMP(k,L(trial,:)==k)).^2) ./ sum(GFP(trial,:).^2);
    end
end

%% 对整个trail的样本按照1s的数据计算4个时间统计参数
list = size(L,2);
if list == 47001
    Lcut = L(:,1:47000);
    Xcut = X(:,1:47000);
elseif list == 46601
    Lcut = L(:,1:46600);
    Xcut = X(:,1:46600);
elseif list == 41201
    Lcut = L(:,1:41200);
    Xcut = X(:,1:41200);
elseif list == 47601
    Lcut = L(:,1:47600);
    Xcut = X(:,1:47600);
elseif list == 37001
    Lcut = L(:,1:37000);
    Xcut = X(:,1:37000);
elseif list == 39001
    Lcut = L(:,1:39000);
    Xcut = X(:,1:39000);
elseif list == 47401
    Lcut = L(:,1:47400);
    Xcut = X(:,1:47400);
elseif list == 43201
    Lcut = L(:,1:43200);
    Xcut = X(:,1:43200);
elseif list ==53001
    Lcut = L(:,1:53000);
    Xcut = X(:,1:53000);
end
m = size(Lcut,2);
s = 1;
for i = 1:200:m
    sampleL = Lcut(:,i:i+199);
    sampleX = Xcut(:,i:i+199);
%% Microstate order for transition probabilities 微状态序列的跃迁概率
order = cell(Ntrials,1);
for trial = 1:Ntrials
    [order{trial}, ~ ] = my_RLE(sampleL(trial,:));
end


%% Preallocating arrays for stats and readying GFP
% prepare arrays (only MGFP can have NANs!)
MGFP = nan(Ntrials,K);  % 生成1*5的矩阵
MDur = zeros(Ntrials,K);
MOcc = zeros(Ntrials,K);
TCov = zeros(Ntrials,K);

GFPcut = squeeze(std(sampleX));  % std(X)函数求解的是最常见的标准差,此时除以的是N-1，默认为按列
if Ntrials>1
    GFPcut = GFPcut';
end

%% For each MS class...
for k = 1:K
    for trial = 1:Ntrials
        
        % Mean GFP
        MGFP(trial,k) = nanmean(GFPcut(trial,sampleL(trial,:)==k));  % 函数对矩阵GFP每一列求均值，计算之前忽略A中的NaN元素
        
        [runvalue, runs] = my_RLE(sampleL(trial,:));
        
        % Mean Duration
        if isnan(nanmean(runs(runvalue == k)))
            MDur(trial,k) = 0;
            MOcc(trial,k) = 0;
        else
            MDur(trial,k) =  nanmean(runs(runvalue == k)) .* (1000 / fs);
            % Occurence
            MOcc(trial,k) =  length(runs(runvalue == k))./200.* fs;
        end
        % time coverage
        TCov(trial,k) = (MDur(trial,k) .* MOcc(trial,k))./ 1000;
        
       
    end
end
MGFPsample(s,:) = MGFP;
MDursample(s,:) = MDur;
MOccsample(s,:) = MOcc;
TCovsample(s,:) = TCov;
%% Transition Probabilities (as with hmmestimate(states,states);)
for trial = 1:Ntrials
    states = order{trial,:};
    states = states(states ~= 0);
    % prepare output matrix
    numStates = K;
    tr = zeros(numStates);
    seqLen = length(states);
    % count up the transitions from the state path
    for count = 1:seqLen-1
        tr(states(count),states(count+1)) = tr(states(count),states(count+1)) + 1;
    end
    trRowSum = sum(tr,2);
    % if we don't have any values then report zeros instead of NaNs.
    trRowSum(trRowSum == 0) = -inf;
    % normalize to give frequency estimate.
    TP(:,:,trial) = tr./repmat(trRowSum,1,numStates);
end
TPsample{s,:} = TP;
s = s+1;
end

%% Write to EEG structure
%   per trial:
Mstats.GEVtotal = GEVtotal;
Mstats.Gfp = MGFPsample;
Mstats.Occurence = MOccsample;
Mstats.Duration = MDursample;
Mstats.Coverage = TCovsample;
Mstats.GEV = GEV;
Mstats.MspatCorr = MspatCorr;
Mstats.TP = TPsample;
Mstats.seq = seq;
Mstats.msFirstTF = msFirstTF;
Mstats.polarity = polarity;

if Ntrials > 1
    % mean parameters
    Mstats.avgs.Gfp = nanmean(MGFP,1);
    Mstats.avgs.Occurence = nanmean(MOcc,1);
    Mstats.avgs.Duration = nanmean(MDur,1);
    Mstats.avgs.Coverage = nanmean(TCov,1);
    Mstats.avgs.GEV = nanmean(GEV,1);
    Mstats.avgs.MspatCorr = nanmean(MspatCorr,1);
    
    % standard deviation of parameters 参数标准差
    Mstats.avgs.stdGfp = nanstd(MGFP);
    Mstats.avgs.stdOccurence = nanstd(MOcc);
    Mstats.avgs.stdDuration = nanstd(MDur);
    Mstats.avgs.stdCoverage = nanstd(TCov);
    Mstats.avgs.stdGEV = nanstd(GEV);
    Mstats.avgs.stdMspatCorr = nanstd(MspatCorr);
end
end

function [d,c]=my_RLE(x)
%% RLE
% This function performs Run Length Encoding to a strem of data x.该函数对数据流x执行运行长度编码
% [d,c]=rl_enc(x) returns the element values in d and their number of  返回d中的元素值及其出现次数放入c中，
% apperance in c. All number formats are accepted for the elements of x.  x中的元素可以接受所有数字格式。
% This function is built by Abdulrahman Ikram Siddiq in Oct-1st-2011 5:15pm.

if nargin~=1
    error('A single 1-D stream must be used as an input') % 必须使用单个一维流作为输入
end

ind=1;
d(ind)=x(1);
c(ind)=1;

for i=2 :length(x)
    if x(i-1)==x(i)
        c(ind)=c(ind)+1;
    else ind=ind+1;
        d(ind)=x(i);
        c(ind)=1;
    end
end

end
