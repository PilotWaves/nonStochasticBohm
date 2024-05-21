tic

sample_size = 100000;
speed = 0.1;
wall = 5;
wall_width = speed;
boundary = 15;

set(0,'DefaultAxesFontSize',20)

    x = randn(1,sample_size);
    x = (x-min(x))/(max(x) - min(x)) - (boundary-wall);
    v=0*x;

    omega = 95;
    nu = 1/2.1;
    K = omega/2/pi;
    C = 1/3;

    a = rand(1,sample_size);
    b = rand(1,sample_size);

    writerObj = VideoWriter('Tunneling.avi');
    writerObj.FrameRate = 10; 
    open(writerObj);

while max(x) < wall
    
    b = 4*b.*(1-b);
    a = mod(a + b*sqrt(2),1);

    v = C*(v + K*cos(2*pi*a).*sin(x).*exp(-nu*abs(v)));
    x = x + v/50 + speed;

    if max(x) > -10
    histogram(x, 'Normalization','pdf','edgecolor','none')
    hold on
    plot([wall+speed/2 wall+speed/2],[0 2.5],'--','linewidth',1)
    hold off
%     pd = fitdist(x', 'Normal');
%     hold on
%     plot([min(x):0.001:max(x)],pdf(pd,[min(x):0.001:max(x)]),'linewidth',4)
%     hold off
    axis([-10 15 0 2.5])
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    end
end

count = 0;
reflect = 0*x;
past = reflect;

while max(x) < boundary

    count = count + 1;

    for i = 1:length(x)
        if reflect(i) == 1
            b(count+1,i) = 4*b(count,i)*(1-b(count,i));
            a(count+1,i) = mod(a(count,i) + b(count+1,i)*sqrt(2),1);

            v(count+1,i) = C*(v(count,i) + K*cos(2*pi*a(count+1,i))*sin(x(count,i))*exp(-nu*abs(v(count,i))));
            x(count+1,i) = x(count,i) - v(count+1,i)/50  - speed;

        else
            b(count+1,i) = 4*b(count,i)*(1-b(count,i));
            a(count+1,i) = mod(a(count,i) + b(count+1,i)*sqrt(2),1);

            v(count+1,i) = C*(v(count,i) + K*cos(2*pi*a(count+1,i))*sin(x(count,i))*exp(-nu*abs(v(count,i))));
            x(count+1,i) = x(count,i) + v(count+1,i)/50 + speed;

            if x(count+1,i) >= wall && x(count+1,i) < wall+wall_width
                x(count+1,i) = 2*wall - x(count+1,i);
                reflect(i) = 1;
            end

            if x(count+1,i) >= wall+wall_width
                past(i) = 1;
            end

        end
    end

    if reflect + past > 0*reflect
        
        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot([wall+speed/2 wall+speed/2],[0 2.5],'--','linewidth',1)
        hold off
        axis([-10 15 0 2.5])
        frame = getframe(gcf);
        writeVideo(writerObj,frame);

    else

        histogram(x(count+1,:),'Normalization','pdf','edgecolor','none')
        hold on
        plot([wall+speed/2 wall+speed/2],[0 2.5],'--','linewidth',1)
        hold off
        axis([-10 15 0 2.5])
        for repeat = 1:10
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end

    end

end

close(writerObj)

toc