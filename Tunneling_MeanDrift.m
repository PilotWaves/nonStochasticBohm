tic

sample_size = 100000;
speed = 0.1;
wall = 5;
wall_height = speed;
wall_width = 1/2;
boundary = 13;

set(0,'DefaultAxesFontSize',20)

    x = randn(1,sample_size);
    x = (x-min(x))/(max(x) - min(x)) - (boundary-wall);
    v=0*x;

    nu = 1/2.1;
    K = 10;
    C = 1/3;

    a = rand(1,sample_size);
    b = rand(1,sample_size);

    writerObj = VideoWriter('Tunneling_MeanDrift.avi');
    writerObj.FrameRate = 10; 
    open(writerObj);

while max(x) < wall
    
    b = 4*b.*(1-b);
    a = mod(a + b*sqrt(2),1);

    v = C*(v + K*cos(2*pi*a).*sin(x).*exp(-nu*abs(v)));
    x = x + v/50 + speed;
    [kde, xi] = ksdensity(x);

    tiledlayout(2,1,'TileSpacing','compact','Padding','tight')

    if max(x) > 0
    
    nexttile
    histogram(x, 'Normalization','pdf','edgecolor','none')
    hold on
    patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], 'red', 'FaceAlpha',0.2)
    hold on
    plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
    hold off
%     pd = fitdist(x', 'Normal');
%     hold on
%     plot([min(x):0.001:max(x)],pdf(pd,[min(x):0.001:max(x)]),'linewidth',4)
%     hold off
    axis([-5 15 0 2])
    xticks([])
    

    nexttile
    histogram(x, 'Normalization','pdf','edgecolor','none')
    hold on
    plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
    hold off
    axis([-5 15 0 2])

    frame = getframe(gcf);
    writeVideo(writerObj,frame);
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

    tiledlayout(2,1,'TileSpacing','compact','Padding','tight')

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

        nexttile    
        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], 'red', 'FaceAlpha',0.2)
        hold on
        plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold on
        plot([xfit.mean xfit.mean], [0 max(pdf(xfit,[min(x_past):0.001:max(x_past)]))], '--', 'linewidth', 2)
        a1 = annotation('textbox',[xfit.mean/15 0.75 0.1 0.1],'String',xfit.mean,'FitBoxToText','on', 'EdgeColor', 'none');
        hold off
        axis([-5 15 0 2])
        xticks([])

        nexttile
        histogram(x_sans(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot(xi_sans,kde_sans,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold on
        plot([xfit_sans.mean xfit_sans.mean], [0 max(pdf(xfit_sans,[min(x_sans(count+1,:)):0.001:max(x_sans(count+1,:))]))], '--', 'linewidth', 2)
        a2 = annotation('textbox',[xfit_sans.mean/15 0.25 0.1 0.1],'String',xfit_sans.mean,'FitBoxToText','on', 'EdgeColor', 'none');
        hold off
        axis([-5 15 0 2])

        frame = getframe(gcf);
        writeVideo(writerObj,frame);

        delete(a1)
        delete(a2)

    else

        nexttile
        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        patch([wall wall+wall_width wall+wall_width wall],[0 0 wall_height/speed wall_height/speed], 'red', 'FaceAlpha',0.2)
        hold off
        axis([-5 15 0 2])
        xticks([])

        nexttile
        histogram(x_sans(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot(xi_sans,kde_sans,'Color',[0 0.4470 0.7410],'linewidth',4)
        hold off
        axis([-5 15 0 2])

        for repeat = 1:10
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end

    end

end

close(writerObj)

toc