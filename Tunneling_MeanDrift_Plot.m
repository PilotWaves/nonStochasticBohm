tic

sample_size = 100000;
speed = 0.1;
wall = 5;
wall_height = speed;
wall_width = 1/2;
boundary = 15;

set(0,'DefaultAxesFontSize',20)

    x = randn(1,sample_size);
    x = (x-min(x))/(max(x) - min(x)) - (boundary-wall);
    v=0*x;

    nu = 1/2.1;
    K = 10;
    C = 1/3;

    a = rand(1,sample_size);
    b = rand(1,sample_size);

    frame_count = 0;
    tiledlayout(3,2,"TileSpacing","compact",'Padding','tight')

while max(x) < wall
    
    b = 4*b.*(1-b);
    a = mod(a + b*sqrt(2),1);

    v = C*(v + K*cos(2*pi*a).*sin(x).*exp(-nu*abs(v)));
    x = x + v/50 + speed;
    [kde, xi] = ksdensity(x);

    if max(x) > 0 

        if frame_count == 0
    
    nexttile
    histogram(x, 'Normalization','pdf','edgecolor','none')
    hold on
    plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
    hold on
    patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], [0.8500 0.3250 0.0980], 'FaceAlpha',0.5)
    hold off
%     pd = fitdist(x', 'Normal');
%     hold on
%     plot([min(x):0.001:max(x)],pdf(pd,[min(x):0.001:max(x)]),'linewidth',4)
%     hold off
    axis([-3 13 0 2])
    xticks([])
    

    nexttile
    histogram(x, 'Normalization','pdf','edgecolor','none')
    hold on
    plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
    hold off
    axis([-3 13 0 2])
    xticks([])
    yticks([])

        end

    frame_count = frame_count + 1;

    end
end

count = 0;
reflect = 0*x;
past = reflect;
v_sans = v;
x_sans = x;

while max(x) < boundary

    count = count + 1;

    for i = 1:sample_size
        if reflect(i) == 1
            b(count+1,i) = 4*b(count,i)*(1-b(count,i));
            a(count+1,i) = mod(a(count,i) + b(count+1,i)*sqrt(2),1);

            v(count+1,i) = C*(v(count,i) + K*cos(2*pi*a(count+1,i))*sin(x(count,i))*exp(-nu*abs(v(count,i))));
            x(count+1,i) = x(count,i) - v(count+1,i)/50  - speed;

            v_sans(count+1,i) = C*(v_sans(count,i) + K*cos(2*pi*a(count+1,i))*sin(x_sans(count,i))*exp(-nu*abs(v_sans(count,i))));
            x_sans(count+1,i) = x_sans(count,i) + v_sans(count+1,i)/50 + speed;

        else
            b(count+1,i) = 4*b(count,i)*(1-b(count,i));
            a(count+1,i) = mod(a(count,i) + b(count+1,i)*sqrt(2),1);

            v(count+1,i) = C*(v(count,i) + K*cos(2*pi*a(count+1,i))*sin(x(count,i))*exp(-nu*abs(v(count,i))));
            x(count+1,i) = x(count,i) + v(count+1,i)/50 + speed;

            v_sans(count+1,i) = C*(v_sans(count,i) + K*cos(2*pi*a(count+1,i))*sin(x_sans(count,i))*exp(-nu*abs(v_sans(count,i))));
            x_sans(count+1,i) = x_sans(count,i) + v_sans(count+1,i)/50 + speed;

            if x(count+1,i) >= wall && x(count+1,i) < wall+wall_height
                x(count+1,i) = 2*wall - x(count+1,i);
                reflect(i) = 1;
            end

            if x(count+1,i) >= wall+wall_width
                past(i) = 1;
            end

        end
    end

            [kde, xi] = ksdensity(x(count+1,:));
            [kde_sans, xi_sans] = ksdensity(x_sans(count+1,:));
            xfit_sans = fitdist(x_sans(count+1,:)','Normal');


    if reflect + past > 0*reflect

        x_past = 0;
        count_past = 0;

        for i = 1:sample_size
            if past(i) == 1
                count_past = count_past + 1;
                x_past(count_past) = x(count+1,i);
            end
        end

        xfit = fitdist(x_past','Normal');

        if frame_count == 114
        nexttile    
        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold on
        plot([xfit.mean xfit.mean], [0 max(pdf(xfit,[min(x_past):0.001:max(x_past)]))], '--', 'linewidth', 2)
        a1 = annotation('textbox',[0.35 0.15 0.1 0.1],'String',round(xfit.mean, 4, 'significant'),'FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 18);
        hold on
        patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], [0.8500 0.3250 0.0980], 'FaceAlpha',0.5)
        hold off
        axis([-3 13 0 2])

        nexttile
        histogram(x_sans(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot(xi_sans,kde_sans,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold on
        plot([xfit_sans.mean xfit_sans.mean], [0 max(pdf(xfit_sans,[min(x_sans(count+1,:)):0.001:max(x_sans(count+1,:))]))], '--', 'linewidth', 2)
        a2 = annotation('textbox',[0.83 0.15 0.1 0.1],'String',round(xfit_sans.mean, 4, 'significant'),'FitBoxToText','on', 'EdgeColor', 'none', 'FontSize', 18);
        hold off
        axis([-3 13 0 2])
        yticks([])
        %delete(a1)
        %delete(a2)
        frame_count
        end

        frame_count = frame_count + 1;

    else

        if frame_count == 64
        nexttile
        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], [0.8500 0.3250 0.0980], 'FaceAlpha',0.5)
        hold off
        axis([-3 13 0 2])
        xticks([])

        nexttile
        histogram(x_sans(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot(xi_sans,kde_sans,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold off
        axis([-3 13 0 2])
        xticks([])
        yticks([])
        %pause
        %frame_count
        end

        frame_count = frame_count + 1;

    end

end

%exportgraphics(gcf,'Tunneling_MeanDrift.pdf')

toc