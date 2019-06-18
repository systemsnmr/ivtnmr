function save_figure(fig_handle,file_name)

set(fig_handle,'Units','Inches');
pos = get(fig_handle,'Position');
set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

if ~exist('figure_images', 'dir')
  mkdir('figure_images');
end

print(fig_handle,fullfile('figure_images',file_name),'-dpdf','-r1200')

% CONVENIENT TO USE with mfilename variable - which returns the name of the current file!
% 
% -r option sets the resolution (also for the VECTOR IMAGES!) -- set to 300-1200 if contour approximation is bad.
% -r0
% -r300
% -r1200
end