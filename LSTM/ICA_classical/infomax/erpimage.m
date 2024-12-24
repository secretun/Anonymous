% erpimage() - image single-trial ERPs optionally sorted on and/or aligned to 
%                 an input variable and smoothed by moving-average (Note: to
%                 return event-aligned data without plotting, use eventlock()).
%                 Click on axes to examine separately and zoom.
% Usage:
%   >> [outdata,outvar,outtrials,limits,axhndls,...
%         erp,amps,cohers,cohsig,ampsig,...
%           outamps,phsangls,phsamp,sortidx] = ...
%             erpimage(data,sortvar,times,'title',avewidth,decimate,[option(s)]);
% Inputs:
%   data     - single-channel data: format (1,frames*trials) or (frames,trials)
%   sortvar  - vector variable to sort trials on (ntrials = length(sortvar))
%                     Example: sortvar = rts (in ms)
%   times    - vector of times in ms (frames=length(times)){def|0->[0:frames-1]}
%              OR [startms ntimes srate] = start time (ms), sampling rate (Hz),
%                                          ntimes = time points per epoch
%  'title'   - title string {default none}
%   avewidth - ntrials in moving average window (may be non-int) {def|0->1}
%   decimate - factor to decimate ntrials out by (may be non-int) {def|0->1}
%               If this is large (> sqrt(number of trials)), this many trials output.
%   option(s)- 'align',[time] -> time lock data to sortvar aligned to time in msec
%                 (time = Inf -> align to median sortvar) {default: no align}
%            - 'nosort'-> don't sort data on sortvar {default: sort}
%            - 'noplot'-> don't plot sortvar {default: plot if in times range}
%            - 'limits',[lotime hitime minerp maxerp loamp hiamp locoher hicoher] 
%                        (can use nan for missing items and omit late items)
%            - 'caxis',[lo hi] -> set color axis limits {default: data bounds}
%                   or [fraction] = set caxis limits at (+/-)fraction*max(abs(data))
%            - 'cbar' -> plot color bar to right of erp-image {default no}
%            - 'erp' -> plot erp time average of the trials below the image
%            - 'phase', [time prctile freq] -> sort by phase at freq (Hz)
%                        in window ending at given time, where
%                        prctile = [0-100] percent data to reject for low amp.
%                        or [-100-0] high amplitude {default: [0 25 8 13]}
%                    or [time prct minfrq maxfrq] -> sort by phase at max power
%                        frequency in the data within the range [minfrq,maxfrq] 
%                    or [time prct minfrq maxfrq topphase] -> sort by phase 
%                        and put topphase (-180<=deg<=180) at top {default 180}
%            - 'coher',[freq] -> plot erp plus amp & coher at freq (Hz)
%                    (or at phase freq above, if specified).
%                or [freq alpha] -> add coher. signif. level line at alpha
%                   (alpha range: (0,0.1]) {default none}
%            - 'allamps' -> image amplitudes at each time & trial. Requires arg
%                   'coher' with alpha signif. {default: image raw data}
%            - 'allcohers',[data2] -> image the coherences at each time & trial. 
%                   Requires arg 'coher' with alpha significance. 
%                   Shows projection on grand mean coherence vector at each time 
%                   and trial. {default: no}
%            - 'topo',{map,eloc_file} -> plot a 2-d head map (vector) at upper left. 
%                   See >> topoplot('example') for electrode location file structure.
%            - 'spec',[loHz,hiHz] -> plot the mean data spectrum at upper right. 
%            - 'srate',[freq]-> specify data sampling rate in Hz 
%                         (if not given in times arg above)
%            - 'vert',[times] -> plot vertical dotted lines at specified times
%            - 'noxlabel' -> do not plot "Time (ms)" on the bottom x-axis
% Outputs:
%   outdata  = (times,epochsout) data matrix (after smoothing)
%   outvar   = (1,epochsout)  sortvar vector (after smoothing)
%   outtrials= (1,epochsout)  smoothed trial numbers 
%   limits   = (1,9) array, as in option 'limits' above plus analysis frequency
%   axhndls  = vector of 1-4 plot axes handles
%   erp      = plotted ERP average
%   amps     = mean amplitude time course
%   coher    = mean inter-trial phase coherence time course
%   cohsig   = coherence significance level
%   ampsig   = amplitude significance levels [lo high]
%   outamps  = matrix of imaged amplitudes
%   phsangls = vector of sorted trial phases at the phase-sorting frequency
%   phsamp   = vector of sorted trial amplitudes at the phase-sorting frequency
%   sortidx  = indices of sorted data epochs plotted

% Tzyy-Ping Jung & Scott Makeig, CNL / Salk Institute, La Jolla 3-2-98
% Uses external toolbox functions: phasecoher(), rmbase(), cbar()
% Uses included functions:         plot1erp(), phasedet(),
%
% 3/5/98 added nosort option -sm
% 3/22/98 added colorbar ylabel, sym. range finding -sm
% 5/08/98 added noplot option -sm
% 6/09/98 added align, erp, coher options -sm
% 6/10/98 added limits options -sm
% 6/26/98 made number of variables output 8, as above -sm 
% 9/16/98 plot out-of-bounds sortvars at nearest times boundary -sm
% 10/27/98 added cohsig, alpha -sm
% 10/28/98 adjust maxamp, maxcoh computation -sm
% 05/03/99 added horizontal ticks beneath coher trace, fixed vert. 
%          scale printing -t-pj
% 05/07/99 converted amps plot to log scaling -sm
% 05/14/99 added sort by alpha phase -se
% 07/23/99 made "Phase-sorted" axis label -sm
% 07/24/99 added 'allamps' option -sm
% 08/04/99 added new times spec., 'srate' arg, made 'phase' and 'allamps'
%          work together, plot re-aligned time zeros  -sm
% 06/26/99 debugged cbar; added vert lines at aligntime to plot1erp() axes -sm
% 09/29/99 fixed srate computation from times -sm & se
% 01/18/00 output outsort without clipping -sm
% 02/29/00 added 'vert' arg, fixed xticklabels, added ampsig -sm
% 03/03/00 added 'vert' arg lines to erp/amp/coher axes -sm
% 03/17/00 added axcopy -sm & tpj
% 03/29/00 added 'topo' option -sm 
% 05/05/00 fixed y-axis label bug when time limits given -sm
% 06/01/00 added topphase arg to 'phase' option for phasemovie.m -sm
% 07/12/00 adjusted prctle()
% 07/12/00 added 'spec' option -sm
% 08/22/00 added coherfreq to limits output -sm
% 09/13/00 added hard limit (1) to maxcoh -sm
% 09/14/00 made topoplot() and psd() plots relative to gca, not gcf -sm
% 10/10/00 added NoTimeflag -sm
% 11/03/00 changed method of rejecting small amplitude trials for phase sort -sm
% 11/03/00 added number_of_trials_out option for decfactor -sm
% 11/16/00 added ampoffset to center sig lines around baseline mean amp (0) -sm
%
% Known Bugs:
% 'limits', [lotime hitime] does not work with 'erp'
% 'limits', [... loerp hierp] (still??) may duplicate "ghost" grey numbers on the coher axis?

function [data,outsort,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,allamps,phaseangles,phsamp,sortidx] = erpimage(data,sortvar,times,titl,avewidth,decfactor,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)

%
%%%%%%%%%%%%%%%%%%%%%%% Define defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
YES = 1;  % logical variables
NO  = 0;

TIMEX = 1;          % 1 -> plot time on x-axis; 0 -> trials on x-axis
BACKCOLOR = [0.8 0.8 0.8]; % default grey background
SORTWIDTH = 1.3;    % 2;
ZEROWIDTH = 1.6;    % 1.6
VERTWIDTH = 3.0;    % 2.0
PLOT_HEIGHT = 0.2;  % fraction of y dim taken up by each time series axes
YGAP = 0.03;        % fraction gap between time axes
YEXPAND = 1.3;      % expansion factor for y-axis about erp, amp data limits
DEFAULT_AVEWIDTH =  1; % smooth trials with this window size by default
DEFAULT_DECFACTOR = 1; % decimate by this factor by default
DEFAULT_CYCLES = 3; % use this many cycles in amp,coher computation window
DEFAULT_CBAR = NO;  % do not plot color bar by default
DEFAULT_PHARGS = [0 25 8 13]; % Default arguments for phase sorting
DEFAULT_ALPHA = 0.01;

LABELFONT = 14;     % font sizes for axis labels, tick labels
TICKFONT  = 11;
alpha     = 0;      % default alpha level for coherence significance

Noshow    = NO;     % show sortvar by default
Nosort    = NO;     % sort on sortvar by default
Caxflag   = NO;     % use default caxis by default
Caxis     = [];
caxfraction = [];
Coherflag = NO;     % don't compute or show amp,coher by default
Cohsigflag= NO;     % default: do not compute coherence significance
Allampsflag=NO;     % don't image the amplitudes by default
Allcohersflag=NO;   % don't image the coherence amplitudes by default
Topoflag  = NO;     % don't plot a topoplot in upper left
Specflag  = NO;     % don't plot a spectrum in upper right
Erpflag   = NO;     % don't show erp average by default
Alignflag = NO;     % don't align data to sortvar by default
Colorbar  = NO;     % if YES, plot a colorbar to right of erp image
Limitflag = NO;     % plot whole times range by default
Phaseflag = NO;     % don't sort by alpha phase
Srateflag = NO;     % srate not given
Vertflag  = NO;
verttimes = [];
NoTimeflag = NO;    % by default DO print "Time (ms)" below bottom axis
coherfreq = nan;    % amp/coher-calculating frequency
freq = 0;           % phase-sorting frequency
srate     = nan;
aligntime = nan;
timelimits= nan;
topomap   = [];     % topo map vector
lospecHz  = [];     % spec lo frequency
topphase = 180;     % default top phase for 'phase' option

minerp = nan; % default limits
maxerp = nan;
minamp = nan;
maxamp = nan;
mincoh = nan;
maxcoh = nan;
baseamp =nan;

ax1    = nan; % default axes handles
axcb   = nan;
ax2    = nan;
ax3    = nan;
ax4    = nan;

%
%%%%%%%%%%%%%%%%%%% Test, fill in commandline args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin < 2
  help erpimage
  return
end

framestot = size(data,1)*size(data,2);
ntrials = length(sortvar);
if ntrials < 2
  help erpimage
  fprintf('\nerpimage(): too few trials.\n');
  return
end

frames = floor(framestot/ntrials);
if frames*ntrials ~= framestot
  help erpimage
  fprintf(...
    '\nerpimage(); length of sortvar doesnt divide no. of data elements.\n')
  return
end

if nargin < 6
  decfactor = 0;
end
if nargin < 5
  avewidth = 0;
end
if nargin<4
  titl = '';
end
if nargin<3
  times = NO;
end
if length(times) == 1 | times == NO,
   times = 0:frames-1;
   srate = 1000*(length(times)-1)/(times(length(times))-times(1));
   fprintf('Using sampling rate %g Hz.\n',srate);
elseif length(times) == 3
   mintime = times(1);
   frames = times(2);
   srate = times(3);
   times = mintime:1000/srate:mintime+(frames-1)*1000/srate;
   fprintf('Using sampling rate %g Hz.\n',srate);
else
   srate = 1000*(length(times)-1)/(times(end)-times(1));
end
if length(times) ~= frames
   fprintf(...
'erpimage(): length(data)(%d) ~= length(sortvar)(%d) * length(times)(%d).\n\n',...
                  framestot,              length(sortvar),   length(times));
   return
end
if avewidth == 0,
  avewidth = DEFAULT_AVEWIDTH;
end
if decfactor == 0,
  decfactor = DEFAULT_DECFACTOR;
end
if avewidth < 1
  help erpimage
  fprintf('\nerpimage(): Variable avewidth cannot be < 1.\n')
  return
end
if avewidth > ntrials
  fprintf('Setting variable avewidth to max %d.\n',ntrials)
  avewidth = ntrials;  
end
if decfactor < 1
  help erpimage
  fprintf('\nerpimage(): Variable decfactor cannot be < 1.\n')
  return
end
if decfactor > ntrials
  fprintf('Setting variable decfactor to max %d.\n',ntrials)
  decfactor = ntrials;  
end

if decfactor > sqrt(ntrials) % if large, output this many trials
  n = 1:ntrials;
  if exist('phargs') & length(phargs)>1
    if phargs(2)>0
      n = n(ceil(phargs(2)*ntrials)+1:end); % trials after rejection
    elseif phargs(2)<0
      n = n(1:floor(phargs(2)*length(n)));  % trials after rejection
    end
  end
  decfactor = length(n)/decfactor;
end
%
%%%%%%%%%%%%%%%%% Collect option args %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin > 6
  flagargs = [];

  for a=7:nargin % for each remaining Arg

    Arg = eval(['arg' int2str(a-6)]);
    if Caxflag == YES
      if size(Arg,1) ~= 1 | size(Arg,2) > 2
        help erpimage
        fprintf('\nerpimage(): caxis arg must be a scalar or (1,2) vector.\n');
        return
      end
      if size(Arg,2) == 2
         Caxis = Arg;
      else
         caxfraction = Arg;
      end
      Caxflag = NO;
    elseif Coherflag == YES
       if size(Arg,1) ~= 1 | ( size(Arg,2) ~= 1 & size(Arg,2) ~= 2)
        help erpimage
        fprintf('\nerpimage(): coher arg must be size (1,1) or (1,2).\n');
        return
       end
       coherfreq = Arg(1);
       if size(Arg,2) == 2
         Cohsigflag = YES;
         alpha  = Arg(2);
         if alpha < 0 | alpha > 0.1
           fprintf('erpimage(): alpha value %g out of bounds.\n',alpha); 
           return
         end
       end
       Coherflag = NO;
       Erpflag = YES;  % plot amp, coher below erp time series
    elseif Topoflag == YES;
       if length(Arg) ~= 2
        help erpimage
        fprintf('\nerpimage(): topo arg must be a list of length 2.\n');
        return
       end
       topomap = Arg{1};
       eloc_file = Arg{2};
       if ~exist(eloc_file)
        fprintf('\nerpimage(): electrode locations file not found.\n');
        return
       end
       Topoflag = NO;
    elseif Specflag == YES;
       if length(Arg) ~= 2
        help erpimage
        fprintf('\nerpimage(): spec arg must be a list of length 2.\n');
        return
       end
       lospecHz = Arg(1);
       hispecHz = Arg(2);
       Specflag = NO;
    elseif Alignflag == YES
       if length(Arg) ~= 1 
        help erpimage
        fprintf('\nerpimage(): align arg must be a scalar msec.\n');
        return
       end
       aligntime = Arg(1);
       Alignflag = NO;
    elseif Limitflag == YES
      %  [lotime hitime loerp hierp loamp hiamp locoher hicoher]
      if size(Arg,1) ~= 1 | size(Arg,2) < 2 ...
          | size(Arg,2) > 9 ...
        help erpimage
        fprintf('\nerpimage(): limits arg must be a vector sized (1,2<->9).\n');
        return
      end
      if  ~isnan(Arg(1)) & (Arg(2) <= Arg(1))
        help erpimage
        fprintf('\nerpimage(): time limits out of order or out of range.\n');
        return
      end
      if Arg(1) < min(times)
          Arg(1) = min(times);
          fprintf('Adjusting mintime limit to first data value %g\n',min(times));
      end
      if Arg(2) > max(times)
          Arg(2) = max(times);
          fprintf('Adjusting maxtime limit to last data value %g\n',max(times));
      end
      timelimits = Arg(1:2);
      if length(Arg)> 2
        minerp = Arg(3);
      end
      if length(Arg)> 3
        maxerp = Arg(4);
      end
      if ~isnan(maxerp) & maxerp <= minerp
        help erpimage
        fprintf('\nerpimage(): erp limits args out of order.\n');
        return
      end
      if length(Arg)> 4
        minamp = Arg(5);
      end
      if length(Arg)> 5
        maxamp = Arg(6);
      end
      if maxamp <= minamp
        help erpimage
        fprintf('\nerpimage(): amp limits args out of order.\n');
        return
      end
      if length(Arg)> 6
        mincoh = Arg(7);
      end
      if length(Arg)> 7
        maxcoh = Arg(8);
      end
      if maxcoh <= mincoh
        help erpimage
        fprintf('\nerpimage(): coh limits args out of order.\n');
        return
      end
      if length(Arg)>8
        baseamp = Arg(9);
      end
      Limitflag = NO;

     elseif Srateflag == YES
          srate = Arg(1);
          Srateflag = NO;
     elseif Vertflag == YES
          verttimes = Arg;
          Vertflag = NO;
     elseif Allcohersflag == YES
          data2=Arg;
          if size(data2) ~= size(data)
fprintf('erpimage(): allcohers data matrix must be the same size as data.\n');
              return
          end
          Allcohersflag = NO;
     elseif Phaseflag == YES
          n = length(Arg);
          if n > 5
            error('erpimage(): Too many arguments for keyword ''phase''');
          end
          phargs = Arg;

          if phargs(3) < 0
              error('erpimage(): Invalid negative argument for keyword ''phase''');
          end
          if n>=4
            if phargs(4) < 0
              error('erpimage(): Invalid negative argument for keyword ''phase''');
            end
          end
          
          if min(phargs(1)) < times(1) | max(phargs(1)) > times(end)
            error('erpimage(): End time for phase sorting filter out of bound.');
          end

          if phargs(2) >= 100 | phargs(2) < -100
            error('%-argument for keyword ''phase'' must be (-100;100)');
          end
          
          if length(phargs) >= 4 & phargs(3) > phargs(4)
            error('erpimage(): Phase sorting frequency range must be increasing.');
          end
          if length(phargs) == 5
            topphase = phargs(5);
          end
          Phaseflag = NO;
    elseif strcmp(Arg,'nosort')
      Nosort = YES;
    elseif strcmp(Arg,'noplot')
      Noshow = YES;
    elseif strcmp(Arg,'caxis')
      Caxflag = YES;
    elseif strcmp(Arg,'coher')
      Coherflag = YES;
    elseif strcmp(Arg,'allamps')
      Allampsflag = YES;
    elseif strcmp(Arg,'allcohers')
      Allcohersflag = YES;
    elseif strcmp(Arg,'topo') | strcmp(Arg,'topoplot')
      Topoflag = YES;
    elseif strcmp(Arg,'spec') | strcmp(Arg,'spectrum')
      Specflag = YES;
    elseif strcmp(Arg,'erp')| strcmp(Arg,'ERP')
      Erpflag = YES;
    elseif strcmp(Arg,'align')
      Alignflag = YES;
    elseif strcmp(Arg,'cbar') | strcmp(Arg,'colorbar')
      Colorbar = YES;
    elseif strcmp(Arg,'limits')
      Limitflag = YES;
    elseif strcmp(Arg,'phase')
      Phaseflag = YES;
    elseif strcmp(Arg,'srate')
      Srateflag = YES;
    elseif strcmp(Arg,'vert')
      Vertflag = YES;
    elseif strcmp(Arg,'noxlabel') | strcmp(Arg,'noxlabels') | strcmp(Arg,'nox')
      NoTimeflag = YES;
    else
      help erpimage
      if isstr(Arg)
         fprintf('\nerpimage(): unknown arg %s\n',Arg);
      else
         fprintf('\nerpimage(): unknown arg %d, size(%d,%d)\n',a,size(Arg,1),size(Arg,2));
      end
      return
    end
  end % Arg
end

if   Caxflag == YES ...
  |Coherflag == YES ...
  |Alignflag == YES ...
  |Limitflag == YES
    help erpimage
    fprintf('\nerpimage(): missing option arg.\n')
    return
end
if (Allampsflag | exist('data2')) & ( isnan(coherfreq) | ~Cohsigflag )
 fprintf('\nerpimage(): allamps and allcohers flags require coher freq, srate, and cohsig.\n');
 return
end
if Allampsflag & exist('data2')
 fprintf('\nerpimage(): cannot image both allamps and allcohers.\n');
 return
end
if ~exist('srate') | srate <= 0 
  fprintf('\nerpimage(): Data srate must be specified and > 0.\n');
  return
end
if exist('phargs')
 if phargs(3) > srate/2
  fprintf('erpimage(): Phase-sorting frequency must be less than Nyquist rate.');
 end
 if length(phargs)==4 & phargs(4) > srate/2
    phargs(4) = srate/2;
 end
 if length(phargs)==5 & (phargs(5)>180 | phargs(5) < -180)
  fprintf('\nerpimage(): coher topphase (%g) out of range.\n',topphase);
  return
 end
end
if ~isnan(coherfreq)
 if coherfreq <= 0 | coherfreq > srate/2 | srate <= 0
  fprintf('\nerpimage(): coher frequency (%g) out of range.\n',coherfreq);
  return
 end
end
          
if isnan(timelimits)
   timelimits = [min(times) max(times)];
end
if ~isnan(aligntime)
 if ~isinf(aligntime) ...
      & (aligntime < timelimits(1) | aligntime > timelimits(2))
  help erpimage
  fprintf('\nerpimage(): requested align time outside of time limits.\n');
  return
 end
end
% 
%%%%%%%%%%%%%%%%  Replace nan's with 0s %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
nans= find(isnan(data));
if length(nans)
  fprintf('Replaced %d nan in data with 0s.\n');
  data(nans) = 0;
end
%
%%%%%%%%%%%%%% Reshape data to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if size(data,2) ~= ntrials
   if size(data,1)>1
% fprintf('frames %d, ntrials %d length(data) %d\n',frames,ntrials,length(data));
     data=reshape(data,1,frames*ntrials);
   end
   data=reshape(data,frames,ntrials);
end
fprintf('Plotting input data as %d epochs of %d frames sampled at %3.1f Hz.\n',...
                             ntrials,frames,srate);
%
%%%%%%%%%%%%%% Reshape data2 to (frames,ntrials) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('data2')
  if size(data2,2) ~= ntrials
   if size(data2,1)>1
     data2=reshape(data2,1,frames*ntrials);
   end
   data2=reshape(data2,frames,ntrials);
  end
end
%
%%%%%%%%%%%%%%%%%%% Align data to sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(aligntime)
  if isinf(aligntime)
    ssv = sort(sortvar); % ssv = 'sorted sortvar'
    aligntime= median(sortvar);
      % Alternative: trimmed median - ignore top/bottom 5%
      % aligntime= median(ssv(ceil(ntrials/20)):floor(19*ntrials/20)); 
    fprintf('Aligning data to median sortvar.\n'); 
  end
  fprintf('Realigned sortvar plotted at %g ms.\n',aligntime);

  aligndata=0*ones(frames,ntrials); % begin with matrix of nan's
  shifts = zeros(1,ntrials);
  for t=1:ntrials, %%%%%%%%% foreach trial %%%%%%%%%
   shft = round((aligntime-sortvar(t))*srate/1000);
   shifts(t) = shft;
   if shft>0, % shift right
    if frames-shft > 0
     aligndata(shft+1:frames,t)=data(1:frames-shft,t);
    else
     fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
    end
   elseif shft < 0 % shift left
    if frames+shft > 0
     aligndata(1:frames+shft,t)=data(1-shft:frames,t);
    else
     fprintf('No aligned data for epoch %d - shift (%d frames) too large.\n',t,shft);
    end
   else % shft == 0
     aligndata(:,t) = data(:,t);
   end 
  end % end trial
  fprintf('Shifted epochs by %d to %d frames.\n',min(shifts),max(shifts));
  data = aligndata;                       % now data is aligned to sortvar
end 
%
%%%%%%%%%%%%%%% Sort the data on sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if exist('phargs') == 1
        % fprintf('     length(phargs) = %d\n',length(phargs))
        if length(phargs) >= 4 % find max frequency in specified band
                [pxx,freqs] = psd(data(:),1024,srate,frames,0);

           %gf = gcf;
           % figure;plot(freqs,pxx);
           %xx=axis;
           %axis([phargs(3) phargs(4) xx(3) xx(4)]);
           %figure(gf);

           pxx = 10*log10(pxx);
           n = find(freqs >= phargs(3) & freqs <= phargs(4));
           if ~length(n)
                freq = phargs(3);
           end
           [dummy maxx] = max(pxx(n));
           freq = freqs(n(maxx));
        else
           freq = phargs(3); % else use specified frequency
        end

        [dummy minx] = min(abs(times-phargs(1)));
        winlen = floor(3*srate/freq);
        winloc = minx-[winlen:-1:0];
        winloc = winloc(find(winloc>0 & winloc<=frames));

        [phaseangles phsamp] = phasedet(data,frames,srate,winloc,freq);

        fprintf(...
  'Sorting data epochs by phase at %.2f Hz, window ending at %f3. ms.\n',...  
             freq,phargs(1));
        fprintf('Phase is computed using a filter of length %d frames.\n',...
             length(winloc));
        %
        % Reject small (or large) phsamp trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        phargs(2) = phargs(2)/100; % convert from % to fraction
        [tmp n] = sort(phsamp);
        if phargs(2)>=0
          n = n(ceil(phargs(2)*length(n))+1:end);
          fprintf(...
    'Retaining %d epochs with largest power at the analysis frequency,\n',...
                      length(n));
          % n = find(prctle(phsamp,phargs(2)) <= phsamp);
        else % phargs(2) < 0
          phargs(2) = 1+phargs(2); % subtract from end
          n = n(1:floor(phargs(2)*length(n)));
          % n = find(prctle(phsamp,phargs(2)) >= phsamp);
          fprintf(...
    'Retaining %d epochs with smallest power at the analysis frequency,\n',...
                      length(n));
        end
        data = data(:,n);
        %
        % Remove low-amplitude trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        phaseangles = phaseangles(n);
        sortvar = sortvar(n);
        ntrials = length(n);
        %
        % Sort remaining data by phase angle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        topphase = topphase/360*2*pi;
        phaseangles = -phaseangles;
        ip = find(phaseangles>topphase);
        phaseangles(ip) = phaseangles(ip)-2*pi;
        [phaseangles sortidx] = sort(phaseangles);
        data = data(:,sortidx);
        phaseangles = -phaseangles;
        phsamp = phsamp(n(sortidx));
        sortvar = sortvar(sortidx);

        fprintf('Size of data = [%d,%d]\n',size(data,1),size(data,2));
        sortidx = n(sortidx); % return original indices in sorted order

elseif Nosort == YES
  fprintf('Not sorting data on input sortvar.\n');
  sortidx = 1:ntrials;	
else
  fprintf('Sorting data on input sortvar.\n');
  [sortvar,sortidx] = sort(sortvar);
  data = data(:,sortidx);
end
if max(sortvar)<0
   fprintf('Inverting sortvar to make it positive.\n');
   sortvar = -sortvar;
end
% 
%%%%%%%%%%%%%%%%%% Smooth data using moving average %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% if ~isnan(coherfreq)
  urdata = data; % compute amp, coher on unsmoothed data
% end
if ~Allampsflag & ~exist('data2') % if imaging potential,
 if avewidth > 1 | decfactor > 1
  if Nosort == YES
    fprintf('Smoothing the data using a window width of %g epochs ',avewidth);
  else
    fprintf('Smoothing the sorted epochs with a %g-epoch moving window.',...
                       avewidth);
  end
  fprintf('\n');
  fprintf('  and a decimation factor of %g\n',decfactor);
  [data,outtrials] = movav(data,1:ntrials,avewidth,decfactor); 
                                            % Note: using square window
  [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
  fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outtrials));
 else % don't smooth
  outtrials = 1:ntrials;
  outsort = sortvar;
 end

 %
 %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 if ~isempty(Caxis) 
  mindat = Caxis(1);
  maxdat = Caxis(2);
  fprintf('Using the specified caxis range of [%g,%g].\n',...
                                           mindat,maxdat);
 else
  mindat = min(min(data));
  maxdat = max(max(data));
  maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
  mindat = -maxdat;
  if ~isempty(caxfraction)
     adjmax = (1-caxfraction)/2*(maxdat-mindat);
     mindat = mindat+adjmax;
     maxdat = maxdat-adjmax;
     fprintf(...
 'The caxis range will be %g times the sym. abs. data range -> [%g,%g].\n',...
                                  caxfraction,mindat,maxdat);
  else
     fprintf(...
 'The caxis range will be the sym. abs. data range -> [%g,%g].\n',...
                                  mindat,maxdat);
  end
 end
end % if ~Allampsflag & ~exist('data2')
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set time limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isnan(timelimits(1))
    timelimits = [min(times) max(times)];
end
fprintf('Data will be plotted between %g and %g ms.\n',timelimits(1),timelimits(2));
%
%%%%%%%%%%%%% Image the aligned/sorted/smoothed data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(coherfreq)       % if plot three time axes
     image_loy = 3*PLOT_HEIGHT;
elseif Erpflag == YES   % elseif if plot only one time axes
     image_loy = 1*PLOT_HEIGHT;
else                    % else plot erp-image only
     image_loy = 0*PLOT_HEIGHT;
end
gcapos=get(gca,'Position');
delete(gca)
if isempty(topomap)
  image_top = 1;
else
  image_top = 0.9;
end
ax1=axes('Position',...
    [gcapos(1) gcapos(2)+image_loy*gcapos(4) ...
     gcapos(3) (image_top-image_loy)*gcapos(4)]);

ind = isnan(data);    % find nan's in data
[i j]=find(ind==1);
if ~isempty(i)
  data(i,j) = 0;      % plot shifted nan data as 0 (=green)
end

if ~Allampsflag & ~exist('data2') %%%%%%%%%%%%%%%% Plot ERP image %%%%%%%%%%%%%%

 if TIMEX
  imagesc(times,outtrials,data',[mindat,maxdat]);% plot time on x-axis
  set(gca,'Ydir','normal');
  axis([timelimits(1) timelimits(2) ...
       min(outtrials) max(outtrials)]);
 else
  imagesc(outtrials,times,data,[mindat,maxdat]); % plot trials on x-axis
  axis([min(outtrials) max(outtrials)...
       timelimits(1) timelimits(2)]);
 end
 hold on
 drawnow
elseif Allampsflag %%%%%%%%%%%%%%%%% Plot allamps instead of data %%%%%%%%%%%%%%

 if freq > 0 
    coherfreq = freq; % use phase-sort frequency
 end
 if alpha>0
   fprintf('Computing and plotting %g coherence significance level...\n',alpha);
   [amps,cohers,cohsig,ampsig,allamps] = ...
     phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,alpha);
   fprintf('Coherence significance level: %g\n',cohsig);
   fprintf('Amplitude significance levels: [%g %g]\n',ampsig(1),ampsig(2));
 else
   [amps,cohers,cohsig,ampsig,allamps] = ...
     phasecoher(urdata,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
 end
 % fprintf('#1 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
 base = find(times<=0);
 if length(base)<2
     base = 1:floor(length(times)/4); % default first quarter-epoch
 end
 amps = 20*log10(amps); % convert to dB
 ampsig = 20*log10(ampsig); % convert to dB
 if isnan(baseamp)
    [amps,baseamp] = rmbase(amps,length(times),base); % remove baseline
 else
    amps = amps - baseamp;
 end
 ampsig = ampsig-baseamp;
 baseall = mean(mean(allamps(base,:)));
 % fprintf('#2 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
 fprintf('Dividing by the mean baseline amplitude %g\n',baseall);
 allamps = allamps./baseall;
 % fprintf('#3 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));

 if avewidth > 1 | decfactor > 1
   % Note: using square window
   if Nosort == YES
    fprintf(...
       'Smoothing the amplitude epochs using a window width of %g epochs ',...
            avewidth);
   else
    fprintf(...
       'Smoothing the sorted amplitude epochs with a %g-epoch moving window.',...
            avewidth);
   end
   fprintf('\n');
   fprintf('  and a decimation factor of %g\n',decfactor);
   %fprintf('4 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
   [allamps,outtrials] = movav(allamps,1:ntrials,avewidth,decfactor); 
                                            % Note: using square window
   %fprintf('5 Size of allamps = [%d %d]\n',size(allamps,1),size(allamps,2));
   [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
   fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outtrials));
 else
  outtrials = 1:ntrials;
  outsort = sortvar;
 end
 allamps = 20*log10(allamps); % convert to dB
 % allamps = allamps';
 %
 %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 if ~isempty(Caxis) 
   mindat = Caxis(1);
   maxdat = Caxis(2);
  fprintf('Using the specified caxis range of [%g,%g].\n',...
                                           mindat,maxdat);
 else
  mindat = min(min(allamps));
  maxdat = max(max(allamps));
  maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
  mindat = -maxdat;
  if ~isempty(caxfraction)
     adjmax = (1-caxfraction)/2*(maxdat-mindat);
     mindat = mindat+adjmax;
     maxdat = maxdat-adjmax;
     fprintf(...
 'The caxis range will be %g times the sym. abs. data range -> [%g,%g].\n',...
                                  caxfraction,mindat,maxdat);
  else
     fprintf(...
     'The caxis range will be the sym. abs. data range -> [%g,%g].\n',...
                                  mindat,maxdat);
  end
 end
 %
 %%%%%%%%%%%%%%%%%%%%% Image amplitudes at coherfreq %%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 fprintf('Plotting amplitudes at freq %g Hz instead of potentials.\n',coherfreq);
 if TIMEX
  imagesc(times,outtrials,allamps',[mindat,maxdat]);% plot time on x-axis
  set(gca,'Ydir','normal');
  axis([timelimits(1) timelimits(2) ...
       min(outtrials) max(outtrials)]);
 else
  imagesc(outtrials,times,allamps,[mindat,maxdat]); % plot trials on x-axis
  axis([min(outtrials) max(outtrials)...
       timelimits(1) timelimits(2)]);
 end
  drawnow
  hold on

elseif exist('data2') %%%%%% Plot allcohers instead of data %%%%%%%%%%%%%%%%%%%%

 if freq > 0 
    coherfreq = freq; % use phase-sort frequency
 end
 if alpha>0
   fprintf('Computing and plotting %g coherence significance level...\n',alpha);
   [amps,cohers,cohsig,ampsig,allcohers] = ...
     crosscoher(urdata,data2,length(times),srate,coherfreq,DEFAULT_CYCLES,alpha);
   fprintf('Inter-Trial Coherence significance level: %g\n',cohsig);
   fprintf('Amplitude significance levels: [%g %g]\n',ampsig(1),ampsig(2));
 else
   [amps,cohers,cohsig,ampsig,allcohers] = ...
     crosscoher(urdata,data2,length(times),srate,coherfreq,DEFAULT_CYCLES,0);
 end
 if ~exist('allcohers')
     fprintf('erpimage(): allcohers not returned....\n')
     return
 end
 allamps = allcohers; % output variable
 % fprintf('Size allcohers = (%d, %d)\n',size(allcohers,1),size(allcohers,2));
 % fprintf('#1 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
 base = find(times<=0);
 if length(base)<2
     base = 1:floor(length(times)/4); % default first quarter-epoch
 end
 amps = 20*log10(amps); % convert to dB
 ampsig = 20*log10(ampsig); % convert to dB
 if isnan(baseamp)
  [amps,baseamp] = rmbase(amps,length(times),base); % remove baseline
 else
   amps = amps - baseamp;
 end
 ampsig = ampsig-baseamp;
 % fprintf('#2 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));

 if avewidth > 1 | decfactor > 1
     % Note: using square window
   if Nosort == YES
    fprintf(...
      'Smoothing the amplitude epochs using a window width of %g epochs '...
            ,avewidth);
   else
    fprintf(...
       'Smoothing the sorted amplitude epochs with a %g-epoch moving window.'...
            ,avewidth);
   end
   fprintf('\n');
   fprintf('  and a decimation factor of %g\n',decfactor);
   % fprintf('4 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));

   [allcohers,outtrials] = movav(allcohers,1:ntrials,avewidth,decfactor); 
                                            % Note: using square window
   % fprintf('5 Size of allcohers = [%d %d]\n',size(allcohers,1),size(allcohers,2));
   [outsort,outtrials] = movav(sortvar,1:ntrials,avewidth,decfactor); 
   fprintf('Output data will be %d frames by %d smoothed trials.\n',...
                          frames,length(outtrials));
 else
  outtrials = 1:ntrials;
  outsort = sortvar;
 end
 %
 %%%%%%%%%%%%%%%%%%%%%%%%% Find color axis limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 if ~isempty(Caxis) 
   mindat = Caxis(1);
   maxdat = Caxis(2);
  fprintf('Using the specified caxis range of [%g,%g].\n',...
                                           mindat,maxdat);
 else
  mindat = -1;
  maxdat = 1
  fprintf(...
     'The caxis range will be the sym. abs. data range [%g,%g].\n',...
                                                     mindat,maxdat);
 end
 %
 %%%%%%%%%%%%%%%%%%%%% Image coherences at coherfreq %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 fprintf('Plotting coherences at freq %g Hz instead of potentials.\n',coherfreq);
 if TIMEX
  imagesc(times,outtrials,allcohers',[mindat,maxdat]);% plot time on x-axis
  set(gca,'Ydir','normal');
  axis([timelimits(1) timelimits(2) ...
       min(outtrials) max(outtrials)]);
 else
  imagesc(outtrials,times,allcohers,[mindat,maxdat]); % plot trials on x-axis
  axis([min(outtrials) max(outtrials)...
       timelimits(1) timelimits(2)]);
 end
  drawnow
  hold on

end %%%%%%%%%%%%%%%%%%%%%%%%%%% End image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(verttimes)
  fprintf('Plotting vertical lines at times: ');
  for vt = verttimes
     fprintf('%g ',vt);
     if isnan(aligntime) % if non-aligned data
       if TIMEX          % overplot vt on image
         plot([vt vt],[0 max(outtrials)],'k:','Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[vt vt],'k:','Linewidth',VERTWIDTH);
       end
     else                % aligned data
       if TIMEX          % overplot realigned vt on image
         plot(aligntime+vt-outsort,outtrials,'k:','LineWidth',VERTWIDTH); 
       else
         plot(outtrials,aligntime+vt-outsort,'k:','LineWidth',VERTWIDTH); 
       end                                                 
     end
  end
  fprintf('\n');
end

set(gca,'FontSize',TICKFONT)
hold on;
%
%%%%%%%%%%% plot vertical line at 0 or align time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(aligntime) % if trials time-aligned 
 if times(1) <= aligntime & times(frames) >= aligntime
  plot([aligntime aligntime],[0 ntrials],'k','Linewidth',2.0); 
     % plot vertical line at aligntime
 end
else % trials not time-aligned 
 if times(1) <= 0 & times(frames) >= 0
  plot([0 0],[0 ntrials],'k:','Linewidth',2.4); % plot vertical line at time 0
 end
end

if Noshow == NO & ( min(outsort) < timelimits(1) ...
                   |max(outsort) > timelimits(2))
  ur_outsort = outsort; % store the pre-adjusted values
  fprintf('Not all sortvar values within time vector limits: \n')
  fprintf('        outliers will be shown at nearest limit.\n');
  i = find(outsort< timelimits(1));
  outsort(i) = timelimits(1);
  i = find(outsort> timelimits(2));
  outsort(i) = timelimits(2);
end

if TIMEX
 if Nosort == YES
  l=ylabel('Trial Number');
 else
  if exist('phargs')
    l=ylabel('Phase-sorted Trials');
    l=ylabel('Trials');
  else
    l=ylabel('Sorted Trials');
  end
 end
else % if switch x<->y axes
 if Nosort == YES & NoTimeflag==NO
  l=xlabel('Trial Number');
 else
  if exist('phargs')
    l=ylabel('Phase-sorted Trials');
    l=ylabel('Trials');
  elseif NoTimeflag == NO
    l=xlabel('Sorted Trials');
  end
 end
end
set(l,'FontSize',LABELFONT);

t=title(titl);
set(t,'FontSize',LABELFONT);

set(gca,'Box','off');
set(gca,'Fontsize',TICKFONT);
set(gca,'color',BACKCOLOR);
if Erpflag == NO & NoTimeflag == NO
  l=xlabel('Time (ms)');
  set(l,'Fontsize',LABELFONT);
end
%
%%%%%%%%%%%%%%%%%%%% Overplot sortvar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Noshow == YES
  fprintf('Not overplotting sorted sortvar on data.\n');

elseif isnan(aligntime) % plot sortvar on un-aligned data

  if Nosort == NO;
    fprintf('Overplotting sorted sortvar on data.\n');
  end
  hold on; 
  if TIMEX      % overplot sortvar
    plot(outsort,outtrials,'k','LineWidth',SORTWIDTH); 
  else
    plot(outtrials,outsort,'k','LineWidth',SORTWIDTH);
  end                                                 
  drawnow
else % plot re-aligned zeros on sortvar-aligned data
  if Nosort == NO;
    fprintf('Overplotting sorted sortvar on data.\n');
  end
  hold on; 
  if TIMEX      % overplot aligned sortvar on image
    plot([aligntime aligntime],[0 ntrials],'k','LineWidth',SORTWIDTH);
  else
    plot([[0 ntrials],aligntime aligntime],'k','LineWidth',SORTWIDTH);
  end
  fprintf('Overplotting realigned times-zero on data.\n');
  hold on; 

  if TIMEX      % overplot realigned 0-time on image
    plot(0+aligntime-outsort,outtrials,'k','LineWidth',ZEROWIDTH); 
  else
    plot(0+outtrials,aligntime-outsort,'k','LineWidth',ZEROWIDTH); 
  end                                                 
  drawnow
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Plot colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Colorbar == YES
   pos=get(ax1,'Position');
   axcb=axes('Position',...
       [pos(1)+pos(3)+0.02 pos(2) ...
        0.03 pos(4)]);
   cbar(axcb,0,[mindat,maxdat]); % plot colorbar to right of image
    drawnow
   axes(ax1); % reset current axes to the erpimage
end
%
%%%%%%%%%%%%%%%%%%%%%%% Compute ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if Erpflag == YES
 axes(ax1); % reset current axes to the erpimage
 xtick = get(ax1,'Xtick');     % remember x-axis tick locations
 xticklabel = get(ax1,'Xticklabel');     % remember x-axis tick locations
 erp=nanmean(data');           % compute erp average, ignoring nan's
 %
 %%%%%% Plot ERP time series below image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 if isnan(maxerp)
  fac = 10;
  maxerp = 0;
  while maxerp == 0
   maxerp = round(fac*YEXPAND*max(erp))/fac; % minimal decimal places
   fac = 10*fac;
  end
 end
 if isnan(minerp)
  fac = 1;
  minerp = 0;
  while minerp == 0
    minerp = round(fac*YEXPAND*min(erp))/fac; % minimal decimal places
    fac = 10*fac;
  end
 end
 limit = [timelimits(1:2) -max(abs([minerp maxerp])) max(abs([minerp maxerp]))];
          
 if ~isnan(coherfreq)
  set(ax1,'Xticklabel',[]);     % remove tick labels from bottom of image
  ax2=axes('Position',...
     [gcapos(1) gcapos(2)+2/3*image_loy*gcapos(4) ...
      gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);
 else
  ax2=axes('Position',...
     [gcapos(1) gcapos(2) ...
      gcapos(3) image_loy*gcapos(4)]);
 end
 plot1erp(ax2,times,erp,limit); % plot ERP
 if ~isnan(aligntime)
   line([aligntime aligntime],[limit(3:4)*1.1],'Color','k'); % x=median sort value
 end

 set(ax2,'Xtick',xtick);        % use same Xticks as erpimage above
 if ~isnan(coherfreq)
   set(ax2,'Xticklabel',[]);    % remove tick labels from ERP x-axis
 else % bottom axis
   set(ax2,'Xticklabel',xticklabel); % add ticklabels to ERP x-axis
 end

 set(ax2,'Yticklabel',[]);      % remove tick labels from left of image
 set(ax2,'YColor',BACKCOLOR);
 if isnan(coherfreq)            % if no amp and coher plots below . . .
  if TIMEX & NoTimeflag == NO
   l=xlabel('Time (ms)');
   set(l,'FontSize',LABELFONT);
  else
   l=ylabel('Time (ms)');
   set(l,'FontSize',LABELFONT);
  end
 end

 if ~isempty(verttimes)
  for vt = verttimes
     if isnan(aligntime)
       if TIMEX      % overplot vt on ERP
         plot([vt vt],[limit(3:4)],'k:','Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[limit(3:4)],'k:','Linewidth',VERTWIDTH);
       end
     else
       if TIMEX      % overplot realigned vt on ERP
         plot(repmat(median(aligntime+vt-outsort),1,2),[limit(3),limit(4)],...
                                  'k:','LineWidth',VERTWIDTH); 
       else
         plot([limit(3),limit(4)],repmat(median(aligntime+vt-outsort),1,2),...
                                  'k:','LineWidth',VERTWIDTH); 
       end                                                 
     end
   end
end

 ydelta = 1/10*(limit(2)-limit(1)); 
 ytextoffset = limit(1)-1.1*ydelta;
 ynumoffset = limit(1)-0.3*ydelta; 

 t=text(ynumoffset,0.7*limit(3), num2str(limit(3)));
 set(t,'HorizontalAlignment','right','FontSize',TICKFONT)

 t=text(ynumoffset,0.7*limit(4), num2str(limit(4)));
 set(t,'HorizontalAlignment','right','FontSize',TICKFONT)

 ynum = 0.7*(limit(3)+limit(4))/2;
 t=text(ytextoffset,ynum,'\muV','Rotation',90);
 set(t,'HorizontalAlignment','center','FontSize',LABELFONT)

 set(ax2,'Fontsize',TICKFONT);
 set(ax2,'Box','off','color',BACKCOLOR);
  drawnow
else
  erp = [];
end
%
%%%%%%%%%%%%%%%%%%%%% Plot amp, coher time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if ~isnan(coherfreq) 
   if freq > 0 
      coherfreq = freq; % use phase-sort frequency
   end
   %
   %%%%%% Plot amp axis below ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   if ~Allampsflag %%%% don't repeat computation
    if Cohsigflag == NO
     fprintf('Computing and plotting amplitude at %g Hz.\n',coherfreq);
     [amps,cohers] = ...
       phasecoher(urdata,size(times,2),srate,coherfreq,DEFAULT_CYCLES);
    else
     fprintf(...
     'Computing and plotting %g coherence significance level at %g Hz...\n',...
                           alpha,                          coherfreq);
     [amps,cohers,cohsig,ampsig] = ...
       phasecoher(urdata,size(times,2),srate,coherfreq,DEFAULT_CYCLES,alpha);
     fprintf('Coherence significance level: %g\n',cohsig);
    end
    amps = 20*log10(amps); % convert to dB
    fprintf('Data amplitude levels: [%g %g]\n',min(amps),max(amps));
    if alpha>0
       ampsig = 20*log10(ampsig); % convert to dB
       fprintf('Data amplitude significance levels: [%g %g]\n',ampsig(1),ampsig(2));
    end
    if isnan(baseamp)
      base = find(times<=0); % use times<0 as default baseline
      if length(base)<2
         base = 1:floor(length(times)/4); % default first quarter-epoch
      end
       [amps,baseamp] = rmbase(amps,length(times),base); % remove baseline
      fprintf('   Removing baseline amplitude of %d dB for plotting.\n',baseamp);
    else
      fprintf('   Removing baseline amplitude of %d dB for plotting.\n',baseamp);
      amps = amps-baseamp;
    end
    if alpha>0
      ampsig = ampsig-baseamp; % remove baseline
    end
   end % ~Allampsflag

   axis('off') % rm ERP axes axis and labels
   ax3=axes('Position',...
     [gcapos(1) gcapos(2)+1/3*image_loy*gcapos(4) ...
      gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);

   if isnan(maxamp) % if not specified
    fac = 1;
    maxamp = 0;
    while maxamp == 0
     maxamp = floor(YEXPAND*fac*max(amps))/fac; % minimal decimal place
     fac = 10*fac;
    end
    maxamp = maxamp + 10/fac;
   end

   if isnan(minamp) % if not specified
    fac = 1;
    minamp = 0;
    while minamp == 0
     minamp = floor(YEXPAND*fac*max(-amps))/fac; % minimal decimal place
     fac = 10*fac;
    end
    minamp = minamp + 10/fac;
    minamp = -minamp;
   end

ampsig
   ampoffset = mean(amps);
ampoffset
[minamp maxamp]
   if Cohsigflag
     if ampsig(1)-ampoffset>maxamp
       maxamp = ampsig(1)-ampoffset;
     end
     if ampsig(2)-ampoffset< minamp
       minamp = ampsig(2)-ampoffset;
     end
   end
   fprintf('    amps range: [%g,%g]\n',minamp,maxamp);
   plot1erp(ax3,times,amps,[timelimits(1) timelimits(2) minamp maxamp]); % plot AMP

   if ~isnan(aligntime)
     line([aligntime aligntime],[minamp maxamp]*1.1,'Color','k'); 
                                                      % x=median sort value
   end
   set(ax3,'Xtick',xtick);
   set(ax3,'Xticklabel',[]);   % remove tick labels from bottom of image
   set(ax3,'Yticklabel',[]);   % remove tick labels from left of image
   set(ax3,'YColor',BACKCOLOR);
   axis('off');
   set(ax3,'Box','off','color',BACKCOLOR);

   if ~isempty(verttimes)
    for vt = verttimes
     if isnan(aligntime)
       if TIMEX      % overplot vt on amp
         plot([vt vt],[minamp maxamp],'k:','Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[minamp maxamp],'k:','Linewidth',VERTWIDTH);
       end
     else
       if TIMEX      % overplot realigned vt on amp
         plot(repmat(median(aligntime+vt-outsort),1,2),[minamp,maxamp],'k:','LineWidth',VERTWIDTH); 
       else
         plot([minamp,maxamp],repmat(median(aligntime+vt-outsort),1,2),'k:','LineWidth',VERTWIDTH); 
       end                                                 
     end
    end
   end

   if Cohsigflag % plot amplitude significance levels
     hold on
      plot([timelimits(1) timelimits(2)],[ampsig(1) ampsig(1)]-ampoffset,'r','linewidth',2);
      plot([timelimits(1) timelimits(2)],[ampsig(2) ampsig(2)]-ampoffset,'r','linewidth',2);
   end

   t=text(ynumoffset,maxamp, num2str(maxamp,3));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,0, num2str(0));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,minamp, num2str(minamp,3));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ytextoffset,(maxamp+minamp)/2,'dB','Rotation',90);
   set(t,'HorizontalAlignment','center','FontSize',LABELFONT);

   axtmp = axis;
   text(1/13*(axtmp(2)-axtmp(1))+axtmp(1), ...
        11/13*(axtmp(4)-axtmp(3))+axtmp(3), ...
        [num2str(baseamp,4) ' dB']);
    drawnow;
   %
   %%%%%% Make coher axis below amp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %
   ax4=axes('Position',...
     [gcapos(1) gcapos(2) ...
      gcapos(3) (image_loy/3-YGAP)*gcapos(4)]);
   if isnan(maxcoh)
    fac = 1;
    maxcoh = 0;
    while maxcoh == 0
     maxcoh = floor(YEXPAND*fac*max(cohers))/fac; % minimal decimal place
     fac = 10*fac;
    end
    maxcoh = maxcoh + 10/fac;
    if maxcoh>1
       maxcoh=1; % absolute limit
    end
   end
   if isnan(mincoh)
     mincoh = 0;
   end
   coh_handle = plot1erp(ax4,times,cohers,[timelimits mincoh maxcoh]); % plot COHER
   if ~isnan(aligntime)
     line([aligntime aligntime],[[mincoh maxcoh]*1.1],'Color','k'); % x=median sort value
   end
   % set(ax4,'Xticklabel',[]);    % remove tick labels from bottom of image
   set(ax4,'Xtick',xtick);
   set(ax4,'Xticklabel',xticklabel);
   set(ax4,'Yticklabel',[]);    % remove tick labels from left of image
   set(ax4,'YColor',BACKCOLOR);

   if ~isempty(verttimes)
    for vt = verttimes
     if isnan(aligntime)
       if TIMEX      % overplot vt on coher
         plot([vt vt],[mincoh maxcoh],'k:','Linewidth',VERTWIDTH);
       else
         plot([0 max(outtrials)],[mincoh maxcoh],'k:','Linewidth',VERTWIDTH);
       end
     else
       if TIMEX      % overplot realigned vt on coher
         plot(repmat(median(aligntime+vt-outsort),1,2),[mincoh,maxcoh],'k:','LineWidth',VERTWIDTH); 
       else
         plot([mincoh,maxcoh],repmat(median(aligntime+vt-outsort),1,2),'k:','LineWidth',VERTWIDTH); 
       end                                                 
     end
    end
   end

   t=text(ynumoffset,0, num2str(0));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ynumoffset,maxcoh, num2str(maxcoh));
   set(t,'HorizontalAlignment','right','FontSize',TICKFONT);

   t=text(ytextoffset,maxcoh/2,'Coh','Rotation',90);
   set(t,'HorizontalAlignment','center','FontSize',LABELFONT);
    drawnow

   if Cohsigflag % plot coherence significance level
     hold on
     plot([timelimits(1) timelimits(2)],[cohsig cohsig],'r','linewidth',2);
   end

   set(ax4,'Box','off','color',BACKCOLOR);
   set(ax4,'Fontsize',TICKFONT);
   if NoTimeflag==NO
     l=xlabel('Time (ms)');
     set(l,'Fontsize',LABELFONT);
   end
   axtmp = axis;
   text(8/13*(axtmp(2)-axtmp(1))+axtmp(1), ...
        8/13*(axtmp(4)-axtmp(3))+axtmp(3), ...
        [num2str(coherfreq,4) ' Hz']);
else
   amps   = [];    % null outputs unless coherfreq specified
   cohers = [];
end
limits = [timelimits(1:2) minerp maxerp minamp maxamp mincoh maxcoh];
axhndls = [ax1 axcb ax2 ax3 ax4];
if exist('ur_outsort')
   outsort = ur_outsort; % restore outsort clipped values, if any
end
if nargout<1
   data = []; % don't spew out data if no args out and no ;
end

%   
%%%%%%%%%%%%%%% plot a topoplot() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
if (~isempty(topomap)) 
    h(12)=axes('Position',...
    [gcapos(1)+0.10*gcapos(3) gcapos(2)+0.86*gcapos(4),...
        0.20*gcapos(3) 0.14*gcapos(4)]);
    % h(12) = subplot('Position',[.10 .86 .20 .14]); 
    fprintf('Plotting a topo map in upper left.\n');
    topoplot(topomap,eloc_file,'electrodes','off')
    axis('square')
end 
%   
%%%%%%%%%%%%%%% plot a spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
SPECFONT = 10;
if (~isempty(lospecHz)) 
    h(13)=axes('Position',...
        [gcapos(1)+0.82*gcapos(3) ...
         gcapos(2)+0.96*gcapos(4),...
         0.15*gcapos(3)*(0.8/gcapos(3))^0.5 ...
         0.10*gcapos(4)*(0.8/gcapos(4))^0.5]);

    % h(13) = subplot('Position',[.75 .88 .15 .10]); 
    fprintf('Plotting the data spectrum in upper right.\n');
    winlength = frames;
    if winlength > 512
       for k=2:5
          if rem(winlength,k) == 0
            break
          end
       end
       winlength = winlength/k;
    end
        
  % [Pxx, Pxxc, F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P)
    [Pxx,Pxxc,F] = psd(reshape(urdata,1,size(urdata,1)*size(urdata,2)),...
                              512,srate,winlength,0,0.05);
    plot(F,10*log10(Pxx));
    goodfs = find(F>= lospecHz & F <= hispecHz);
    maxgfs = max(10*log10(Pxx(goodfs)));
    mingfs = min(10*log10(Pxx(goodfs)));
    axis('square')
    axis([lospecHz hispecHz mingfs-1 maxgfs+1]);
    set(h(13),'Box','off','color',BACKCOLOR);
    set(h(13),'Fontsize',SPECFONT);
    l=ylabel('dB');
    set(l,'Fontsize',SPECFONT);
    if ~isnan(coherfreq)
      hold on; plot([coherfreq,coherfreq],[mingfs maxgfs],'r');
   end
end 

limits = [limits coherfreq];  % add coherfreq to output limits array
%   
%%%%%%%%%%%%%%% turn on axcopy() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
axcopy; % turn on popup zoom windows

%
%%%%%%%%%%%%%%%%%%% function plot1erp() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [plot_handle] = plot1erp(ax,Time,erp,axlimits)

  LINEWIDTH = 2;
  [plot_handle] = plot(Time,erp,'LineWidth',LINEWIDTH); hold on
  axis([axlimits(1:2) 1.1*axlimits(3:4)])
  l1=line([axlimits(1:2)],[0 0],    'Color','k','linewidth',2.0); % y=zero-line
  l2=line([0 0],[axlimits(3:4)*1.1],'Color','k','linewidth',1.6); % x=zero-line

%
%%%%%%%%%%%%%%%%%%% function phasedet() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phasedet() - function used in erpimage.m
%              Constructs a complex filter at frequency freq
%

function [ang,amp,win] = phasedet(data,frames,srate,nwin,freq)
% Typical values:
%   frames = 768;
%   srate = 256;
%   nwin = [200:300];
%   freq = 10;

data = reshape(data,[frames prod(size(data))/frames]);
win = exp(2i*pi*freq(:)*[1:length(nwin)]/srate);
win = win .* repmat(hanning(length(nwin))',length(freq),1);
resp = win * data(nwin,:);
ang = angle(resp);
amp = abs(resp);
if ~exist('allamps')
   allamps = [];
end
%
%%%%%%%%%%%%%% function prctle() %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function prctl = prctle(data,pc);
[prows pcols] = size(pc);
if prows ~= 1 & pcols ~= 1
    error('pc must be a scalar or a vector.');
end
if any(pc > 100) | any(pc < 0)
    error('pc must be between 0 and 100');
end
[i,j] = size(data);
sortdata = sort(data);
if i==1 | j==1 % if data is a vector
  i = max(i,j); j = 1;
  if i == 1,
    fprintf('  prctle() note: input data is a single scalar!\n')
    y = data*ones(length(pc),1); % if data is scalar, return it
    return;
  end
  sortdata = sortdata(:);
end
pt = [0 100*((1:i)-0.5)./i 100];
sortdata = [min(data); sortdata; max(data)];
prctl = interp1(pt,sortdata,pc);
