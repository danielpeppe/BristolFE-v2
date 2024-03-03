function fn_plot_signal(op, res, steps, exp_data, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Plot Final Results
figure;
FntN='Times New Roman';
FntS = 13;
ax1 = gca;
%Get exp time data for aperture (needed for scaling)
aperture = 1:op.aperture_n_els;
aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
%Iterate over sim results
for i = 1:length(params)
    res_sum_dsps = sum(res{1,i}{1}.dsps); %tmp for readability
    %Scale sim data to exp data
    scale_dsp = max(abs(aperture_dsp_data))/max(abs(res_sum_dsps)); %Scale exp response to match sim response
    translate_time = 1.194e-05; %Start of exp response
    %Plot final results
    plot(steps{1,i}{1}.load.time + translate_time, res_sum_dsps*scale_dsp,'DisplayName',string(params(i)));
    hold on
end

%Plot experimental data on top
hold on
%Plot
plot(exp_data.time, aperture_dsp_data, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','target');
hold off

xlabel('Time (s)')
ylabel('Magnitude (-)')
xlim([0 + translate_time,4e-6 + translate_time])
ylim([min(aperture_dsp_data),max(aperture_dsp_data)])
xcorr = 2; ycorr = 2;
w = 30; h = 20;
osx = 1.8; osy = 1.8;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
legend('Location','south'); legend boxoff
hold off

end

