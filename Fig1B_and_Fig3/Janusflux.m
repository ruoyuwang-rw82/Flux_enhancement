clear all

dslist = linspace(100e-9,100e-6,100);

Jw_dcmd_F = zeros(1, length(dslist));
Jw_dcmd_D = zeros(1, length(dslist));
Jw_vmd = zeros(1, length(dslist));

for i = 1:length(dslist)

    [Jw, Jq, Ts, Pv, B] = model0_janus(150e-6, 0.7, 0.45e-6/2, 0.2, 60, 20, 3e3, 0, dslist(i), 'F'); % LMH
    Jw_dcmd_F(i) = Jw;

    [Jw2, Jq2, Ts2, Pv2, B2] = model0_janus(150e-6, 0.7, 0.45e-6/2, 0.2, 60, 20, 3e3, 0, dslist(i), 'D'); % LMH
    Jw_dcmd_D(i) = Jw2;

    [Jw3, Jq3, Ts3, Pv3, B3] = modelvmd0_janus(150e-6, 0.7, 0.45e-6/2, 60, 10e3, 3e3, 0, dslist(i)); % LMH
    Jw_vmd(i) = Jw3;

end

figure

plot(dslist/1e-6, Jw_dcmd_F,'-','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(dslist/1e-6, Jw_dcmd_D,'--','Color', [253,185,18]/255,'LineWidth',6)
hold on
plot(dslist/1e-6, Jw_vmd,'-','Color', [46,117,182]/255,'LineWidth',6)
hold on

xlim([0.1,100])
ylim([0,50])
pbaspect([1.5 1 1])
set(gca, 'FontSize',30,'linewidth', 3)
set(gca, 'XScale', 'log')
xticks([0.1,1,10,100])
set(gca, 'XTickLabel', arrayfun(@num2str, xticks, 'UniformOutput', false));
set(gca, 'YAxisLocation', 'left', 'TickDir', 'out');