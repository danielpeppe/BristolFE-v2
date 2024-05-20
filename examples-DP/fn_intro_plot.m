function fn_intro_plot()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

aperture_n_els = 8;
load('g4-s8_x-8000_z0_025.mat','exp_data')

%Plot Final Results
fig = figure;
FntN = 'Cambria';
FntS = 18;
%Get exp time data for aperture (needed for scaling)
aperture = 1:aperture_n_els;
aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);

%Plot experimental data on top
%Plot
translate_time = 1.1955e-05; %Start of exp response
plot(exp_data.time - translate_time, aperture_dsp_data, 'Color', 'k','LineStyle','-','DisplayName','Response','LineWidth', 1);
hold off

xlabel('Time [s]')
ylabel('Displacement [mm]')
%CONVERGENCE
% % xlim([2.3e-6,3.25e-06])
xlim([0, 3.5e-6])
ylim([-30,30])
%INTERLAMINAR RESIN PLOT??
% % ylim([-1e-7,1e-7])
xcorr = 2; ycorr = 2;
w = 30; h = 15;
osx = 2.4; osy = 2.4;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
% legend('Location','north'); legend boxoff
set(gca,'XTick',[], 'YTick', [])
hold off

print(fig,'-dpng',['-r','1000'], "intro-plot.png")

end

