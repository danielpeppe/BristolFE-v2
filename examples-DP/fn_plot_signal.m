function fn_plot_signal(op, res, steps, exp_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

params = op.params;

%Plot Final Results
fig = figure;
FntN = 'Cambria';
FntS = 18;
% FntS = 30;
ax = gca;
%Get exp time data for aperture (needed for scaling)
aperture = 1:op.aperture_n_els;
aperture_data = ismember(exp_data.tx, aperture) & ismember(exp_data.rx, aperture);
aperture_dsp_data = sum(exp_data.time_data(:, aperture_data), 2);
%Iterate over sim results
alpha = 1;
for i = 1:length(res)
    res_sum_dsps = sum(res{1,i}{1}.dsps); %tmp for readability
    %Scale sim data to exp data
    scale_sim = 1.8;
    scale_exp = max(abs(aperture_dsp_data))/max(abs(res_sum_dsps)); %Scale exp response to match sim response
    %Plot final results
    if iscell(params)
        params_i = params{i};
        n = length(params_i);
        if n == 2
            for ii = 1:(n - 1)
                params_i = strcat(string(params_i(1)),string(params_i(2)));
            end
        end
        str = params_i;
    else
        str = string(params(i));
        %CONVERGENCE
        % % str = "Elements per wavelength: " + string(params(i));
        % % str = string(params(i));

    end
    %INTERLAMINAR RESIN PLOT
    % % if i == 1
    % %     str = 'Without resin layer';
    % %     plot(steps{1,i}{1}.load.time, res_sum_dsps * scale_sim * op.plot_scale_dsps,'DisplayName',str, 'Color', 'k', 'LineStyle','--', 'LineWidth', 1);
    % % elseif i == 2
    % %     str = 'With resin layer';
    % %     plot(steps{1,i}{1}.load.time, res_sum_dsps * scale_sim  * op.plot_scale_dsps,'DisplayName',str, 'Color', 'k', 'LineStyle','-','LineWidth', 1);
    % % end
    %CONVERGENCE
    % % alpha = alpha - 0.2;
    % % plot(steps{1,i}{1}.load.time, res_sum_dsps * scale_exp * op.plot_scale_dsps,'DisplayName',str, 'Color', [0 0 0]+alpha, 'LineStyle','-', 'LineWidth', 2);
    %VALIDATE WITH EXP
    str = 'Simulated response';
    plot(steps{1,i}{1}.load.time, res_sum_dsps * scale_sim * op.plot_scale_dsps,'DisplayName',str, 'Color', 'k', 'LineStyle','-', 'LineWidth', 1);
    hold on
end

% Plot experimental data on top
hold on
%Plot
translate_time = 1.1955e-05; %Start of exp response
plot(exp_data.time - translate_time, aperture_dsp_data/scale_exp, 'Color', hsv2rgb([.95,1,1]),'LineStyle',':','DisplayName','Experimental response','LineWidth', 1);
hold off

xlabel('Time [s]')
ylabel('Displacement [-]')
%CONVERGENCE
% xlim([2.5e-6,3.25e-06])
xlim([0, 3.5e-6])
% ylim([-30,35])
ylim([-0.7e-7,0.7e-7])
yticklabels({});

% % Adjusting axes properties
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

xcorr = 2; ycorr = 2;
w = 30; h = 15;
osx = 1.1; osy = 2;
% osx = 2; osy = 3.5;
set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
% set(gca,'Units','centimeters','FontName',FntN,'fontsize',FntS)
set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS)
% set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
legend('Location','north'); legend boxoff
% % lgd = legend('Location','northwest'); legend boxoff
% % title(lgd, 'Elements per wavelength:')
hold off

% print(fig,'-dpng',['-r','1000'], "Figures/interlaminar-reflections.png")
% print(fig,'-dpng',['-r','1000'], "Figures/convergence.png")
print(fig,'-dpng',['-r','1000'], "Figures/validate-with-exp.png")

end

