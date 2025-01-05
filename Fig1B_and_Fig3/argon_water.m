dH = 40650;

Pvs = @(x) exp(23.5377-4016.3622./(x-38.6339));
T = 293:1:353;
Pv = Pvs(T);

Tr = 373;
P_pre = Pvs(Tr)*exp(-dH/8.314*(1./T-1/Tr));
P_adj = 1.03e5*exp(dH/8.314/373)*exp(-dH*(1-0.08)/8.314./T);
P_adj2 = 1.03e5*exp(dH/8.314/373)*exp(-dH*(1-0.15)/8.314./T);

figure
plot(T,Pv/1e3,'k--','LineWidth',6)
hold on
plot(T,P_pre/1e3,'-','Color', [46,117,182]/255,'LineWidth',6)
hold on
plot(T,P_adj/1e3,'-','Color', [46,117,182,200]/255,'LineWidth',6)
hold on
plot(T,P_adj2/1e3,'-','Color', [46,117,182,150]/255,'LineWidth',6)
hold on
set(gca, 'YScale', 'log')
pbaspect([1 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');
xlim([280,360])
xticks([280 300 320 340 360])
yticks([1 10 100 1000])
yticklabels({'1', '10', '100', '1000'}) % Specify the corresponding labels


% argon
data_argon = [142.59, 350.59
163.76, 149.10
185.61, 82.03];

H_argon = data_argon(:,1);
density_argon = data_argon(:,2);

Hspan_argon = linspace(5e3,8e3,100);
den_argon_model = 1e5/8.314/90*exp(75/8.314)*exp(-Hspan_argon/8.314/90);

% water
data_water = [[1733.3826363636365, 1945.6716933333335, 2023.876966666667, 2100.2825866666667, 2194.669936170213]', [1.3733095317170376, 0.7138077050813683, 0.3627667303221521, 0.24998451511125874, 0.17563580205911655]'];

H_water = data_water(:,1)*1000*0.018; % J/mol
density_water = data_water(:,2)/0.018; % mol/m3

Hspan_water = linspace(3e4,4.5e4,100);
den_water_model = density_water(end)*exp(H_water(end)/8.314/353)*exp(-Hspan_water/8.314/353);


figure
% Argon data
plot(Hspan_argon/1e3, den_argon_model * 8.314 * 90/1e3, 'b--', ...
    'Color', [253, 185, 18] / 255, 'LineWidth', 6);
hold on
scatter(H_argon * 40 / 1e3, density_argon * 8.314 * 90/1e3, 35^2, ...
    'MarkerEdgeColor', 'k', ...                % Black edge color
    'MarkerFaceColor', [253, 185, 18] / 255, ... % Custom fill color (yellow-orange)
    'LineWidth', 2, ...                        % Line width for the edge
    'MarkerFaceAlpha', 0.75);                   % Transparency for the fill
hold on

% Water data
plot(Hspan_water / 1e3, den_water_model * 8.314 * 353/1e3, 'b--', ...
    'Color', [46, 117, 182] / 255, 'LineWidth', 6);
hold on
scatter(H_water / 1e3, density_water * 8.314 * 353/1e3, 35^2, ...
    'MarkerEdgeColor', 'k', ...                % Black edge color
    'MarkerFaceColor', [46, 117, 182] / 255, ... % Custom fill color (blue)
    'LineWidth', 2, ...                        % Line width for the edge
    'MarkerFaceAlpha', 0.75);                   % Transparency for the fill
hold on
pbaspect([1 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([5,50])
ylim([1,1e3])
yticks([1 10 100 1000])
yticklabels({'1', '10', '100', '1000'}) % Specify the corresponding labels
xticks([5,10,20,30,40,50])
set(gca, 'XTickLabel', arrayfun(@num2str, xticks, 'UniformOutput', false));
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');
