tic

sample_size = 100000;
speed = 0.1;

set(0,'DefaultAxesFontSize',20)
tiledlayout(2,1, "TileSpacing","compact", "Padding","tight")
nexttile

    x = randn(1,sample_size);
    x = x/(max(x) - min(x));
    min_x = min(x);
    v=0*x;
    r = 4;

    nu = 1/2.1;
    K = 10;
    C = 1/3;

    a = rand(1,sample_size);
    b = rand(1,sample_size);

%     writerObj = VideoWriter('FreeSpace.avi');
%     writerObj.FrameRate = 10; 
%     open(writerObj);

hold on
    histogram(x, 'Normalization','pdf','edgecolor','none')
    pd = fitdist(x', 'Normal')
    plot([min(x):0.001:max(x)],pdf(pd,[min(x):0.001:max(x)]),'linewidth',2, 'Color',[0 0.4470 0.7410])
    tickpi = [0, 4, 10];
    set(gca,'XTick',tickpi)
    ticktxt = {'0', '5', '20'};
    set(gca,'XTickLabel',ticktxt)
    %str = {'\Delta x = ' + string(round(max(x)-min(x),2,'significant'))};
    %annotation('textbox', [.15, .7, 1, 1], 'String', str, 'fontsize', 18, 'EdgeColor','none','VerticalAlignment', 'bottom');

    while min(x)<20

    b = r*b.*(1-b);
    a = mod(a + b*sqrt(2),1);

    v = C*(v + K*cos(2*pi*a).*sin(x).*exp(-nu*abs(v)));
    x = x + v/50 + speed;


    if mean(x) < 5.05 && mean(x) > 4.95
    histogram(x-1, 'Normalization','pdf','edgecolor','none')
    pd = fitdist(x'-1, 'Normal')
    plot([min(x-1):0.001:max(x-1)],pdf(pd,[min(x-1):0.001:max(x-1)]),'linewidth',2, 'Color', [0.9290 0.6940 0.1250])
    %str = {'\Delta x = ' + string(round(max(x)-min(x),2,'significant'))};
    %annotation('textbox', [.37, .45, 1, 1], 'String', str, 'fontsize', 18, 'EdgeColor','none','VerticalAlignment', 'bottom');
    %frame = getframe(gcf);
    %writeVideo(writerObj,frame);
    end

if mean(x) < 20.05 && mean(x) > 19.95
    break
end

    end

    histogram(x-10, 'Normalization','pdf','edgecolor','none')
    pd = fitdist(x'-10, 'Normal')
    plot([min(x-10):0.001:max(x-10)],pdf(pd,[min(x-10):0.001:max(x-10)]),'linewidth',2, 'Color', [0.4660 0.6740 0.1880])
    %str = {'\Delta x = ' + string(round(max(x)-min(x),2,'significant'))};
    %annotation('textbox', [.675, .3, 1, 1],  'String', str, 'fontsize', 18, 'EdgeColor','none','VerticalAlignment', 'bottom');
    %annotation('arrow', [0.45, 0.55], [0.6, 0.6], 'linewidth', 4, 'headlength', 20, 'headwidth', 20);
    %annotation('textbox', [.43, .6, 1, 1],  'String', '\Delta v = 0.2', 'fontsize', 18, 'EdgeColor','none','VerticalAlignment', 'bottom')
    max_x = max(x)-10;
    axis([min_x, max_x, 0, 3.6])
    hold off
    xlabel('x')
    ylabel('PDF')
    %exportgraphics(gcf,'FreeSpace.pdf')


%close(writerObj)

toc