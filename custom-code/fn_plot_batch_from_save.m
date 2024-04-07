clear;
close all;
% restoredefaultpath;
addpath("../code");

figure_batch_plot = figure;

FntN='Times New Roman';
FntS = 13;

% Count the number of .mat files
dir_name = 'C:/Users/danjm/Documents/IRP_data/data/batch10';
mat_files = dir(fullfile(dir_name, '*.mat'));
n_mat_files = numel(mat_files);

%Iterate over sim results
for i = 1:n_mat_files
   
    load([dir_name '/response' char(string(i)) '.mat'], 'dsp_data', 'time_data', 'op_config')

    res_norm = dsp_data/max(abs(dsp_data)); %Scale exp response to match sim response
    
    str = strcat(string(op_config.porosity), ", ", string(op_config.attenuation));
    plot(time_data, res_norm, 'DisplayName', str);
    hold on
end


xlabel('Time (s)')
ylabel('Magnitude (-)')
xlim([0, op_config.max_time])
ylim([-1,1]/3)
xcorr = 2; ycorr = 2;
w = 30; h = 15;
osx = 1.8; osy = 1.8;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
% legend('Location','south'); legend boxoff
hold off
