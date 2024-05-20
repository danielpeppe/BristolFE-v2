function fn_plot_porosity_correlations_report(varargin)

if ~isempty(varargin)
    op_save = varargin{1};
else
    batch_number = 100;
    fprintf('op_save is empty...using batch %d\n', batch_number)

    % Count the number of .mat files
    batch_path = "data/batch" + batch_number;
    mat_files = dir(fullfile(batch_path, '*.mat'));
    n_mat_files = numel(mat_files);
    
    %Iterate over sim results
    op_save = cell(1, n_mat_files);
    for i = 1:n_mat_files
        response_path = batch_path + "/response" + i + ".mat";
        load(response_path, 'op_config')
        op_save{1,i} = op_config;
    end
end

figure1 = figure;

FntN='Times New Roman';
FntS = 16;

x = zeros(size(op_save));
y = zeros(size(op_save));
for i = 1:numel(op_save)
    x(i) = op_save{i}.porosity;
    y(i) = 20*log10(op_save{i}.attenuation) / (2*op_save{i}.specimen_size*1000);
    % y(i) = (20*log10(op_save{i}.attenuation)*6 - 89) / (2*op_save{i}.specimen_size*1000);
    % y(i) = op_save{i}.attenuation;
end
[x_sorted, x_i] = sort(x);
y_sorted = y(x_i);


% New code to fit and plot the best-fit line
p = polyfit(x_sorted, y_sorted, 1);  % p(1) is the slope, p(2) is the intercept
y_fit = polyval(p, x_sorted);  % Generate y values based on the fit
plot3 = plot(x_sorted, y_fit, 'k', 'DisplayName', 'Best-fit line','LineStyle','-','LineWidth', 1);  % Plot the best-fit line
% Print out the slope and intercept in the MATLAB command window
fprintf('The slope (gradient) of the best-fit line is: %f\n', p(1));
fprintf('The y-intercept of the best-fit line is: %f\n', p(2));

hold on

% yyaxis left
plot1 = scatter(x_sorted, y_sorted, 'r', 'filled','DisplayName','Simulation');
% ylim([5.85952 8.61])
ylim([1.9 2.35])
set(gca, 'YColor', 'k');
% ylabel('Equivalent attenuation, $\alpha_{eq}$ [-]', 'Interpreter', 'latex')



% yyaxis right
% plot2 = plot([0 3],[0.5 3],'DisplayName','Target relationship', 'Color', 'k', 'LineStyle','-', 'LineWidth', 0.5);
% ylim(get(gca, 'YLim') + 0.5);
% ylim([0.5 3])
% set(gca, 'YColor', 'k');
ylabel('Attenuation coefficient, $\alpha$ [dB/mm]', 'Interpreter', 'latex')
xlabel('Porosity, $p_v$ [\%]', 'Interpreter', 'latex')

% xlim([0, 3])

% Coordinates to highlight
x_highlight = [1.5 3 3 1.5];
y_limits = ylim;
y_highlight = [y_limits(1) y_limits(1) y_limits(2) y_limits(2)];
fill(x_highlight, y_highlight, 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % Red highlight
text(2.25, (y_limits(1) + y_limits(2))/2, 'Extrapolated region', 'HorizontalAlignment', 'center', 'FontSize', 13, 'FontName', FntN, 'FontAngle','italic');

% ylim([5.5, 8.8])
% xcorr = 2; ycorr = 2;
% w = 15; h = 10;
% osx = 2; osy = 2;
% set(gcf,'Units','centimeters','Position',[xcorr ycorr w+1.75*osx h+1.75*osy])
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k',})
% set(gca,'Units','centimeters','Position',[osx osy w h],'FontName',FntN,'fontsize',FntS,'XMinorTick','on','YMinorTick','on')
set(gca,'Units','centimeters','FontName',FntN,'fontsize',FntS)

% legend(); legend boxoff
str1 = 'Simulation';
str3 = 'Best-fit line';
legend([plot1 plot3],str1,str3,'Location','northwest'); legend boxoff
% handlevec = [plot1 plot2];
% str2 = 'Target relationship';
% legend(handlevec,str1,str2,'Location','northwest'); legend boxoff
hold off

print(figure1,'-dpng',['-r','1000'], 'Figures/porosity_vs_attenuation.png')

end