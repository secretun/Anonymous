%  plotdata() - plot concatenated multichannel data epochs in two-column format
%
% Usage:
%       >> plotdata(data)
%       >> plotdata(data,frames)
%       >> plotdata(data,frames,limits,title,channames,colors,rtitle,ydir) 
%
%   data      = data consisting of consecutive epochs of (chans,frames) 
%   frames    = time frames/points per epoch {0 -> data length}
%  [limits]   = [xmin xmax ymin ymax]  (x's in ms) 
%               {0 (or both y's 0) -> use data limits)
%  'title'    = plot title {0 -> none}
%  'channames'= file of channel names, 4 chars per line 
%                ELSE vector of channel numbers (ex: [33:64])
%                (NB: In file, '.' = space)  {0 -> [1:chans]}
%  'colors'   = file of color codes, 3 chars per line  
%                 ('.' = space) {0 -> default color order}
%  'rtitle'   = right-side plot title {0 -> none}
%  ydir       = y-axis polarity (pos-up = 1; neg-up = -1) {def -> pos-up}
%
% Also: See plottopo, timtopo, envtopo, headplot, compmap, eegmovie

% Scott Makeig & Tzyy-Ping Jung, CNL / Salk Institute, La Jolla CA
% 5-1-96 from showerps.m  -sm from showerp.m -tpj
% 5-3-96 added default channel numbering, frames & title -sm
% 5-17-96 added nargin tests below -sm
% 5-21-96 added right titles -sm
% 6-29-96 removed Postscript file query -sm
% 7-22-96 restored lines to fill printed page with figure -sm
% 7-29-96 added [channumbers] option for channames argument. -sm
%         changed "channels" argument to "channames" in help statement above -sm
% 1-6-97  debugged min/max time and +/- printing -tpj & sm
% 3-3-97  restored previous Default axis parameters at end -sm
% 4-2-97  debugged 32-epoch plotting -sm
% 5-10-97 added no-args check -sm
% 5-20-97 read icadefs.m for MAXPLOTDATACHANS and MAXPLOTDATAEPOCHS -sm
% 6-23-97 use returns instead of errorcodes -sm
% 10-31-97 use normalized PaperUnits for US/A4 compatibility, 
%          fixed [xy]{min,max} printing, added clf, adding Clipping off,
%          fixed scaling, added limits tests -sm & ch
% 11-19-97 removed an 'orient' command that caused problems printing -sm
% 07-15-98 added 'ydir' arg, made pos-up the default -sm
% 07-24-98 fixed 'ydir' arg, and pos-up default -sm
% 01-02-99 added warning about frames not dividing data length -sm
% 02-19-99 debugged axis limits -sm

function plotdata(data,frames,limits,plottitle,channels,colors,righttitle,ydr)

if nargin < 1,
    help plotdata
    return
end

%
% Set defaults
%
FONTSIZE = 16;     % font size to use for labels
TICKFONTSIZE=16;   % font size to use for axis labels
DEFAULT_SIGN = 1;  % default to plotting positive-up (1) or negative-up (-1)
ISSPEC = 0;        % default - 0 = not spectral data, 1 = pos-only spectra

axcolor= get(0,'DefaultAxesXcolor'); % find what the default x-axis color is
plotfile = 'plotdata.ps';
ls_plotfile = 'ls -l plotdata.ps';

%
%%%%%%%%%%%%%%%%%%%%%%%%%% Substitute defaults for missing parameters %%%%%
%
SIGN = DEFAULT_SIGN;
if nargin < 8
   ydr = 0;
end
if ydr == -1
   SIGN = -1;
end
if nargin < 7,
    righttitle = 0; 
end;
if nargin < 6
    colors = 0;
end
if nargin < 5
    channels = [1:size(data,1)]; % default channames = 1:chans
end
if nargin < 4
    plottitle = 0; %CJH
end
limitset = 0;
if nargin < 3,
    limits = 0;
elseif length(limits)>1
    limitset = 1;
end
if nargin < 2,
    frames = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  [chans,framestotal]=size(data);             % data size
  if frames <=0,
    frames = framestotal;    % default
    datasets=1;
  elseif frames==1,
    fprintf('plotdata: cannot plot less than 2 frames per trace.\n');
    return
    datasets=1;
  else
    datasets = fix(framestotal/frames);        % number of traces to overplot
    if datasets*frames < framestotal
         fprintf('\nWARNING: %d points at end of data will not be plotted.\n',...
                    framestotal-datasets*frames);
    end
  end;

  icadefs; % read MAXPLOTDATACHANS constant from icadefs.m

  if chans>MAXPLOTDATACHANS,
    fprintf('plotdata: not set up to plot more than %d channels.\n',...
                       MAXPLOTDATACHANS);
    return
  end;
  if datasets>MAXPLOTDATAEPOCHS 
      fprintf('plotdata: not set up to plot more than %d epochs.\n',...
                       MAXPLOTDATAEPOCHS);
    return
  end;
  if datasets<1
      fprintf('plotdata: cannot plot less than 1 epoch!\n');
      return
  end;
%
%%%%%%%%%%%%% Extend the size of the plotting area in the window %%%%%%%%%%%%
%
  curfig = gcf;
  h=figure(curfig);
  set(h,'Color',BACKCOLOR); % set the background color
  set(h,'PaperUnits','normalized'); % use percentages to avoid US/A4 difference
  set(h,'PaperPosition',[0.0235308 0.0272775 0.894169 0.909249]); % equivalent
  % orient portrait
  axis('normal');
%
%%%%%%%%%%%%%%%%%%%% Read the channel names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if channels==0,
    channels = [1:size(data,1)];
  end;
  if isstr(channels) == 0,
    channames = num2str(channels(:));                   %%CJH
  else,
    chid = fopen(channels,'r');
    if chid <3,
        fprintf('plotdata(): cannot open file %s.\n',channels);
        return
    end;
    channames = fscanf(chid,'%s',[4 MAXPLOTDATACHANS]);
    channames = channames';
       [r c] = size(channames);
    for i=1:r
        for j=1:c
            if channames(i,j)=='.',
                channames(i,j)=' ';
            end;
        end;
    end;
  end; % setting channames
%
%%%%%%%%%%%%%%%%%%%%%%%%% Read the color names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if colors ~=0,
    if ~isstr(colors)
       fprintf('plotdata(): color file name must be a string.\n');
       return
    end
    cid = fopen(colors,'r');
    % fprintf('cid = %d\n',cid);
    if cid <3,
        fprintf('plotdata: cannot open file %s.\n',colors);
        return
    end;
    colors = fscanf(cid,'%s',[3 MAXPLOTDATAEPOCHS]);
    colors = colors';
       [r c] = size(colors);
    for i=1:r
        for j=1:c
            if colors(i,j)=='.',
                colors(i,j)=' ';
            end;
        end;
    end;
  else % use default color order (no yellow!)
     colors =['r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  ';'r  ';'b  ';'g  ';'c  ';'m  '];
     colors = [colors; colors];  % make > 64 available
  end;
  for c=1:length(colors)   % make white traces black unless axis color is white
    if colors(c,1)=='w' & axcolor~=[1 1 1]
         colors(c,1)='k';
    end
  end
%
%%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if limits==0,      % == 0 or [0 0 0 0]
    xmin=0;
    xmax=frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
    yrange = ymax-ymin;
    ymin = ymin - 0.00*yrange;
    ymax = ymax + 0.00*yrange;
  else
    if length(limits)~=4,
      fprintf( ...
       'plotdata: limits should be 0 or an array [xmin xmax ymin ymax].\n');
      return
    end;
    xmin = limits(1);
    xmax = limits(2);
    ymin = limits(3);
    ymax = limits(4);
  end;

  if xmax == 0 & xmin == 0,
    x = (0:1:frames-1);
    xmin = 0;
    xmax = frames-1;
  else
    dx = (xmax-xmin)/(frames-1);
    x=xmin*ones(1,frames)+dx*(0:frames-1); % compute x-values
    xmax = xmax*frames/frames;
  end;
  if xmax<=xmin,
      fprintf('plotdata() - xmax must be > xmin.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
      yrange = ymax-ymin;
      ymin = ymin - 0.00*yrange;
      ymax = ymax + 0.00*yrange;
  end
  if ymax<=ymin,
      fprintf('plotdata() - ymax must be > ymin.\n')
      return
  end

  xlabel = 'Time (ms)';
  if ymin >= 0 & xmin >= 0,          % For all-positive (spectral) data
    ISSPEC = 1;
    SIGN = 1;
    fprintf('\nPlotting positive up. Assuming data are spectra.\n');
    xlabel = 'Freq (Hz)';
    ymin = 0;                        % plot positive-up
  end;

%
%%%%%%%%%%%%%%%%%%%%%%%% Set up plotting environment %%%%%%%%%%%%%%%%%%%%%%%%%
%
  h = gcf;
  % set(h,'YLim',[ymin ymax]);       % set default plotting parameters
  % set(h,'XLim',[xmin xmax]);
  % set(h,'FontSize',18);
  % set(h,'DefaultLineLineWidth',1); % for thinner postscript lines
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Print plot info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  clf;   % clear the current figure

  % print plottitle over (left) subplot 1
  if plottitle==0,
    plottitle = '';
  end
  if righttitle==0,
    righttitle = '';
  end
  subplot(ceil(chans/2),2,1), h=gca;title([plottitle],...
              'FontSize',FONTSIZE); % title plot and
  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);
  set(h,'FontSize',FONTSIZE);            % choose font size

  subplot(ceil(chans/2),2,2), h=gca;title([righttitle],...
              'FontSize',FONTSIZE); % title plot and
  set(h,'FontSize',FONTSIZE);            % choose font size
  set(h,'YLim',[ymin ymax]);            % set default plotting parameters
  set(h,'XLim',[xmin xmax]);

  msg = ['\nPlotting %d traces of %d frames with colors: '];
  for c=1:datasets
    msg = [msg  colors(c,:)];
  end
  msg = [msg ' -> \n'];    % print starting info on screen . . .
  fprintf(...
    '\n  limits: [xmin,xmax,ymin,ymax] = [%4.1f %4.1f %4.2f %4.2f]\n',...
                xmin,xmax,ymin,ymax);
  fprintf(msg,datasets,frames);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

  xdiff=xmax-xmin;
  ydiff=ymax-ymin;

  for P=0:datasets-1, %  for each data epoch
      fprintf('\ntrace %d: ',P+1);

    for I=1:chans,        % for each data channel

      index=(2*((rem(I-1,ceil(chans/2))+1)))-1+floor(2*(I-1)/chans); 
      subplot(ceil(chans/2),2,index); h=gca;    % = 1 3 5 .. 2 4 6 ..
      hold on;                      % plot down left side of page first
      set(h,'YLim',[ymin ymax]);    % set default plotting parameters
      set(h,'XLim',[xmin xmax]);

      axislcolor = get(gca,'Xcolor');   %%CJH

      if P==0 
        if ~ISSPEC
         axis('off');
         plot([0 0],[ymin ymax],'color',axislcolor); % draw vert %%CJH
                                                     % axis at time 0  
        else  % ISSPEC
            axis('off');plot([xmin xmin],[0 ymax],'color',axislcolor); 
        end  
                                                
      % secondx = 200;                               % draw second vert axis 
      % axis('off');plot([secondx secondx],[ymin ymax],'color',axislcolor); 
 
       axis('off');
       plot([xmin xmax],[0 0],'color',axislcolor);   % draw horizontal axis 

       %%%%%%%%%%%%%%%%%%%%%%% Print channel names %%%%%%%%%%%%%%%%%%%%%%%%%%

       if channels~=0,                               % print channames
        if ~ISSPEC
          if ymin <= 0 & ymax >= 0,
                yht = 0;
          else
                yht = mean(SIGN*data(I,1+P*frames:1+P*frames+frames-1));
          end
          axis('off'),h=text(xmin-0.04*xdiff,yht,[channames(I,:)]); 
            set(h,'HorizontalAlignment','right');      % print before traces
            set(h,'FontSize',FONTSIZE);                % choose font size

        % axis('off'),h=text(xmax+0.10*xdiff,yht,[channames(I,:)]);
        %    set(h,'HorizontalAlignment','left');      % print after traces

        else % ISSPEC
          axis('off'),h=text(xmin-0.04*xdiff,ymax/2,[channames(I,:)]); 
            set(h,'HorizontalAlignment','right');      % print before traces
            set(h,'FontSize',FONTSIZE);                % choose font size

        % axis('off'),h=text(xmax+0.10*xdiff,ymax/2,[channames(I,:)]);
        %    set(h,'HorizontalAlignment','left');      % print after traces

        end;
       end; 
      end; 
      %
      %%%%%%%%%%%%%%%%%%%%% Plot two-sided time-series data %%%%%%%%%%%%%%%%%%%
      %
      if ~ISSPEC
                                                    
        plot(x,SIGN*data(I,1+P*frames:1+P*frames+frames-1),colors(P+1));   

        if SIGN > 0
            axis([xmin xmax ymin ymax]);           % set axis bounds (pos up)
        else
            axis([xmin xmax -1*ymax -1*ymin]);     % set axis bounds (neg up)
        end
                                        
        if P==datasets-1,            % on last traces
         if I==floor((chans+1)/2),   % draw +/0 on lowest left plot
            signx = xmin-0.04*xdiff;

          if SIGN > 0  % pos up
            axis('off');hl=text(signx,ymin,num2str(ymin,3));        % text ymin
            axis('off');hi=text(signx,ymax,['+' num2str(ymax,3)]);  % text +ymax
          else         % neg up
            axis('off');hl=text(signx,-1*ymin,num2str(ymin,3));        % text ymin
            axis('off');hi=text(signx,-1*ymax,['+' num2str(ymax,3)]);  % text +ymax
          end
          set(hl,'FontSize',TICKFONTSIZE);         % choose font size
          set(hl,'HorizontalAlignment','right','Clipping','off');
          set(hi,'FontSize',TICKFONTSIZE);         % choose font size
          set(hi,'HorizontalAlignment','right','Clipping','off');
         end

         if I==chans & limitset,    % draw timescale on lowest right plot
            ytick = -ymax-0.25*ydiff;

            tick = [int2str(xmin)]; h=text(xmin,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center',...
                      'Clipping','off');  % center text

            tick = [xlabel]; h=text(xmin+xdiff/2,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center',...
                      'Clipping','off');  % center text

            tick = [int2str(xmax)]; h=text(xmax,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE);         % choose font size
            set(h,'HorizontalAlignment','center',...
                      'Clipping','off');  % center text
         end;
        end;
      %
      %%%%%%%%%%%%%%%%%%%%% Plot positive-up [0,ymax] %%%%%%%%%%%%%%%%%%%%%%%%
      %
      else % ISSPEC
        ymin=0;
        plot(x,SIGN*data(I,1+P*frames:1+P*frames+frames-1),colors(P+1));   
        ymaxm = ymax;
        % ymin = 0.01;
        % ymaxm = 10.^ceil(log(ymax)/log(10.));
        % if ymaxm/2. > ymax,
        %    ymaxm = ymaxm/2.;
        % end;
        axis([xmin xmax ymin ymaxm]);      % set axis values

        if P==datasets-1,                  % on last trace
         if I==floor((chans+1)/2),         % draw +/0 on lowest left plot
          signx = xmin-0.04*xdiff;

          axis('off');h=text(signx,ymax,['+' num2str(ymax,3)]); 
            set(h,'FontSize',TICKFONTSIZE);
            set(h,'HorizontalAlignment','right','Clipping','off');        

          axis('off');h=text(signx,0,'0'); 
            set(h,'FontSize',TICKFONTSIZE);
            set(h,'HorizontalAlignment','right','Clipping','off');    
         end;

         if I==chans,                    % draw freq scale on lowest right plot
            ytick = -0.25*ymax;

          tick = [num2str(round(10*xmin)/10) ]; h=text(xmin,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE);              
            set(h,'HorizontalAlignment','center','Clipping','off'); 

          tick = [xlabel]; h=text(xmin+xdiff/2,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE); 
            set(h,'HorizontalAlignment','center','Clipping','off');

          tick = [num2str(round(10*xmax)/10) ]; h=text(xmax,ytick,tick);
            set(h,'FontSize',TICKFONTSIZE);             
            set(h,'HorizontalAlignment','center','Clipping','off');

        end; % if last chan
      end % if last data
     end; % if ~ISSPEC

     fprintf(' %d',I);
    end; % subplot
  end; % dataset
  fprintf('\n');
%
%%%%%%%%%%%%%%%%%% Make printed figure fill page %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  curfig = gcf;
  h=figure(curfig);
  % set(h,'PaperPosition',[0.2 0.3 7.6 10]); % stretch out the plot on the page
%
%%%%%%%%%%%%%%%%%% Restore plot environment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  % set(h,'DefaultAxesYLim',aylim);      % restore previous plotting parameters
  % set(h,'DefaultAxesXLim',axlim);
  % set(h,'DefaultAxesFontSize',axfont);

if 0,    % START DETOUR XXXXXXXXXXXXX
%
%%%%%%%%%%%%%%%%%% Save plot to disk if asked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   if  plotfile ~= '', %
     n=0; y=1;
     answer = input('plotdata: Save plot as Postscript file? (y/n) ');
     if answer==1,
         fprintf('\nSaving figure as %s ... ',plotfile);
         curfig = gcf;
         h=figure(curfig);
         % set(h,'PaperPosition',[0.2 0.3 7.6 10]); 
                                     % stretch out the plot on the page
         eval (['print -dpsc ' plotfile]);
         fprintf('saved. Move or remove file!\n');
         unix(ls_plotfile);
     end
   end
end       % END DETOUR XXXXXXXXXXXXX
