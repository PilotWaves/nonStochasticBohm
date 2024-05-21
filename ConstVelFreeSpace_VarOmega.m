tic

sample_size = 100000;
speed = 0.1;

    writerObj = VideoWriter('FreeSpace_VarOmega.avi');
    writerObj.FrameRate = 10; 
    open(writerObj);

set(0,'DefaultAxesFontSize',20)


%omega = [1, 100, 1000, 10000];

    x = randn(sample_size, 4);
    x = x./(max(x) - min(x));
    min_x = min(x); min_min_x = min(min_x);
    max_x = max(x);
    pd(1) = fitdist(x(:,1), 'Normal');
    pd(2) = fitdist(x(:,2), 'Normal');
    pd(3) = fitdist(x(:,3), 'Normal');
    pd(4) = fitdist(x(:,4), 'Normal');
    pd_mean = pd(1).mean;
    pd_max(1) = max(pdf(pd(1),[min(x(:,1)):0.001:max(x(:,1))]));
    pd_max(2) = max(pdf(pd(2),[min(x(:,2)):0.001:max(x(:,2))]));
    pd_max(3) = max(pdf(pd(3),[min(x(:,3)):0.001:max(x(:,3))]));
    pd_max(4) = max(pdf(pd(4),[min(x(:,4)):0.001:max(x(:,4))]));
    v=0*x;
    r = 4;
    nu = 1/2.1;
    %K = [0, 15, 15^2, 15^3];
    K = [0, 10, 100, 1000];
    C = 1/3;

    a = rand(sample_size, 4);
    b = rand(sample_size, 4);

    % tiledlayout(2,1, "TileSpacing","compact", "Padding","tight")
    % for i = 2
    %     nexttile
    %     [kde, xi] = ksdensity(x(:,i));
    %         histogram(x(:,i)', 'Normalization','pdf','edgecolor','none')
    %         hold on
    %         plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',2)
    %         axis([-2, 10, 0, pd_max(i)])
    %         hold off
    % 
    %         frame = getframe(gcf);
    %         writeVideo(writerObj,frame);
    % end

    while pd_mean<25

        tiledlayout(2,1, "TileSpacing","compact", "Padding","tight")
        for i = 2
                nexttile

            b(:,i) = r*b(:,i).*(1-b(:,i));
            a(:,i) = mod(a(:,i) + b(:,i)*sqrt(2),1);

            v(:,i) = C*(v(:,i) + K(i)*cos(2*pi*a(:,i)).*sin(x(:,i)).*exp(-nu*abs(v(:,i))));
            x(:,i) = x(:,i) + v(:,i)/50 + speed;

            pd = fitdist(x(:,i), 'Normal');
            pd_mean = pd.mean;

            [kde, xi] = ksdensity(x(:,i));

            histogram(x(:,i)', 'Normalization','pdf','edgecolor','none')
            hold on
            plot(xi,kde,'Color',[0 0.4470 0.7410],'linewidth',2)
            axis([-2, 25, 0, pd_max(i)*1.1])
            hold off
            xlabel('x')
            ylabel('PDF')

        end
    
        frame = getframe(gcf);
        writeVideo(writerObj,frame);

    end

%exportgraphics(gcf,'FreeSpace.pdf')
close(writerObj)

toc