% axcopy() - move, resize, or copy Matlab axes using the mouse
%
% Usage:           >> selaxiscopy
%                  >> selaxiscopy(fig)
%                  >> selaxiscopy('noticks',fig)
%
%     Clicking the left mouse button on an axis and copy
%     the objects in the axis to a new (popup) figure. 
%     Option 'noticks' does not make x and y tickloabelmodes 'auto'
%       in the popup.
%
% 3-16-00 Tzyy-Ping Jung & Scott Makeig, CNL / Salk Institute, La Jolla
% requires copyaxes.m
% 4-2-00 added 'noticks' -sm

function axcopy(fig)

if exist('fig') & strcmp(fig,'noticks')
   noticks = 1;
   if nargin> 1
     shift
   else
     fig = [];
   end
end
if ~exist('fig') | isempty(fig) | fig == 0 
   fig = gcf;
end

hndl= findobj('parent',fig,'type','axes');
offidx=[];
for a=1:length(hndl)                    % make all axes visible
  set(findobj('parent',hndl(a)),'ButtonDownFcn','copyaxis');
end
figure(fig);
set(hndl,'ButtonDownFcn','copyaxis');
%if ~exist('noticks')
  %axis on
  %set(gca,'xticklabelmode','auto')
  %set(gca,'yticklabelmode','auto')
%end
