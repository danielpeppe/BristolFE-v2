function fn_plot_batch(op_save, res, steps, batch_path, save_figure)
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
    print(figure_batch_plot,'-dpng',['-r','1000'], batch_path + "/batch_plot.png")
end

end

