% copyaxis() - helper function for axcopy()
%
%  See: axcopy()

% 3-16-00 Tzyy-Ping Jung, CNL / Salk Institute, La Jolla CA

function copyaxis()

OFF_STD = 0.25; % std dev of figure offset
MIN_OFF = 0.15; % minimum offset for new figure
BORDER  = 0.04;  % screen edge tolerance

fig=gcf;
sel=gca;
%
% Get position for new figure
%
set(fig,'Units','normalized');
place=get(fig,'Position');
axiscolor=get(fig,'Color');
cmap=colormap;
newxy = (OFF_STD*randn(1,2))+place(1,1:2);
newx=newxy(1);newy=newxy(2);
if abs(newx-place(1,1))<MIN_OFF, newx=place(1,1)+sign(newx-place(1,1))*MIN_OFF;end
if abs(newy-place(1,1))<MIN_OFF, newy=place(1,1)+sign(newy-place(1,1))*MIN_OFF;end
if newx<BORDER, newx=BORDER; end
if newy<BORDER, newy=BORDER; end
if newx+place(3)>1-BORDER, newx=1-BORDER-place(3); end
if newy+place(4)>1-BORDER, newy=1-BORDER-place(4); end
newfig=figure('Units','Normalized','Position',[newx,newy,place(1,3:4)]);
%
% Copy object to new figure
%
set(newfig,'Color',axiscolor);
copyobj(sel,newfig);set(gca,'Position',[0.130 0.110 0.775 0.815]);
colormap(cmap);
%
% Increase font size
%
set(findobj('parent',newfig,'type','axes'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',16)
set(get(gca,'YLabel'),'FontSize',16)
set(get(gca,'Title'),'Fontsize',16);
%
% Add xtick and ytick labels if missing
%
if strcmp(get(gca,'Box'),'on')
   set(gca,'xticklabelmode','auto')
   set(gca,'xtickmode','auto')
   set(gca,'yticklabelmode','auto')
   set(gca,'ytickmode','auto')
end
%
% Turn on zoom in the new figure
%
zoom on;
%
% Turn off the axis copy for the new figure
%
hndl= findobj('parent',newfig,'type','axes');
set(hndl,'ButtonDownFcn','','Selected','off');
for a=1:length(hndl) % make all axes visible
  set(findobj('parent',hndl(a)),'ButtonDownFcn','','Selected','off');
end
