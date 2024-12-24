% logimagesc() - make an imagesc(0) plot with log y-axis values (ala semilogy())
%
% Usage:  >> [logfreqs,dataout] = logimagesc(times,freqs,data);
%
% Input:
%         times = vector of x-axis values
%         freqs = vector of y-axis values
%         data = matrix of size (freqs,times)

% Scott Makeig, CNL / Salk Institute, La Jolla CA 4/00
% 08-07-00 made ydir normal -sm

function [lgfreqs,datout] = logimagesc(times,freqs,data)

  if size(data,1) ~= length(freqs)
      fprintf('logfreq(): data matrix must have %d rows!\n',length(freqs));
      datout = data;
      return
  end
  if size(data,2) ~= length(times)
      fprintf('logfreq(): data matrix must have %d columns!\n',length(times));
      datout = data;
      return
  end
  if min(freqs)<= 0
      fprintf('logfreq(): frequencies must be > 0!\n');
      datout = data;
      return
  end
  lfreqs = log(freqs);
  lgfreqs = linspace(lfreqs(1),lfreqs(end),length(lfreqs));
  lgfreqs = lgfreqs(:);
  lfreqs = lfreqs(:);
  [mesht meshf] = meshgrid(times,lfreqs);
  datout = griddata(mesht,meshf,data,times,lgfreqs);

  imagesc(times,freqs,data);
  nt = ceil(min(freqs)); % new tick - round up min freq to int-hz
  yt=get(gca,'ytick');
  yl=get(gca,'yticklabel');

  imagesc(times,lgfreqs,datout);
  set(gca,'ydir','normal')
  set(gca,'ytick',log([nt yt]));
  set(gca,'yticklabel',{int2str(nt) yl});
