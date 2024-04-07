function fn_plot_batch(op_save, res, steps, exp_data, save_figure, data_folder_path)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Take first op_save as op to extract common options
op = op_save{1};

%Plot Final Results
if save_figure
    figure_batch_plot = figure('Visible', 'off');
else
    figure_batch_plot = figure;
end

FntN='Times New Roman';
FntS = 13;

%Iterate over sim results
for i = 1:length(res)
    res_sum = sum(res{1,i}{1}.dsps); %tmp for readability
    time = steps{1,i}{1}.load.time;
    %Scale sim data
    res_norm = res_sum/max(abs(res_sum)); %Scale exp response to match sim response
    %Plot final results
    n_pts_half = round(length(res_norm)/2);
    op_save{1,i}.attenuation = 1/max(res_norm(n_pts_half:end));
    
    str = strcat(string(op_save{1,i}.porosity), ", ", string(op_save{1,i}.attenuation));
    plot(time, res_norm, 'DisplayName', str);
    hold on
end

% %Plot experimental data on top
% hold on
% %Get exp time data for aperture (needed for scaling)
% aperture = 1:op.aperture_n_els;
% aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
% aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
% aperture_dsp_data_norm = aperture_dsp_data/max(abs(aperture_dsp_data));
% translate_time = 1.194e-05; %Start of exp response
% %Plot
% exp_time_translated = exp_data.time - translate_time;
% exp_i = find(exp_time_translated > 0 & exp_time_translated < op.max_time);
% correction = 0.75; %exp displacements need to be reduced slightly because exp initial signal need to be smaller than sim's
% plot(exp_time_translated(exp_i), aperture_dsp_data_norm(exp_i) * correction, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target');
% hold off

xlabel('Time (s)')
ylabel('Magnitude (-)')
xlim([0, op.max_time])
ylim([-1,1]/3)
xcorr = 2; ycorr = 2;
w = 30; h = 15;
osx = 1.8; osy = 1.8;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
legend('Location','south'); legend boxoff
hold off

%% Save figure
if save_figure
    print(figure_batch_plot,'-dpng',['-r','1000'], [data_folder_path '/batch_plot.png'])
end

end

