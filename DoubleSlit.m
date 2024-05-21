tic

sample_size = 100001;
y_hist = 0;
z_hist = 0;
count = 1;
slit_pos = 30;
slit_width = 10;

    N = 1001;
    N_half = (N-1)/2 + 1;

    nu = 1/2.1;
    K = 10;
    C = 1/3;

for i = 1:sample_size

    x = 0:N-1;
    vy = 0*x;
    y = vy;

    vz = vy;

    a = rand;
    b = rand;

    y(1) = 6e-3*i/sample_size - 3e-3;
    z = y;

for n = 1:N_half-1
    
    b = 4*b*(1-b);
    a = mod(a + b*sqrt(2),1);

    vy(n+1) = C*(vy(n) + K*cos(2*pi*a)*sin(y(n))*exp(-nu*abs(vy(n))));
    y(n+1) = y(n) + vy(n+1);
    vz(n+1) = C*(vz(n) + K*cos(2*pi*a)*sin(z(n))*exp(-nu*abs(vz(n))));
    z(n+1) = z(n) + vz(n+1);

end

    if (y(N_half) > -slit_pos && y(N_half) < -slit_pos + slit_width) || (y(N_half) > slit_pos - slit_width && y(N_half) < slit_pos)

%         if y(N_half) > -slit_pos && y(N_half) < -slit_pos + slit_width
%             center = -slit_pos + slit_width/2;
%         end
% 
%         if y(N_half) > slit_pos - slit_width && y(N_half) < slit_pos
%             center = slit_pos - slit_width/2;
%         end

        a = 0;
        b = 0.1;
        
        for n = N_half:N-1

            b = 4*b*(1-b);
            a = mod(a + b*sqrt(2),1);

            vy(n+1) = C*(vy(n) + K*cos(2*pi*a)*sin(y(n))*exp(-nu*abs(vy(n))));
            y(n+1) = y(n) + vy(n+1);
            vz(n+1) = C*(vz(n) + K*cos(2*pi*a)*sin(z(n))*exp(-nu*abs(vz(n))));
            z(n+1) = z(n) + vz(n+1);
        end

    end

    %plot(y, x)
    %hold on

    if ((y(N_half) > -slit_pos && y(N_half) < -slit_pos + slit_width)...
            || (y(N_half) > slit_pos - slit_width && y(N_half) < slit_pos))% ...
            %&& (abs(y(end)) < 5*pi - pi/2)
        y_hist(count) = y(end);
        z_hist(count) = z(end);
        %max(abs(v))
        count = count + 1;
%         plot(x,y, 'linewidth',0.5, 'color', [0.4660 0.6740 0.1880])
%         hold on
%     else
%         plot(x,y,'--', 'linewidth',0.5, 'color', [0.9290 0.6940 0.1250])
%         hold on
    end

end

%     plot([500 500], [-100 -30], 'k', [500 500], [-20 20], 'k', [500 500], [30 100], 'k', 'linewidth',5)
%     hold off
%     axis([0 1000 -100 100])


    %figure(1)
    %plot([-100 -30], [500 500], 'k', [-20 20], [500 500], 'k', [30 100], [500 500], 'k', 'linewidth',4)
    %hold off
    %axis([-100 100 0 1000])



    %figure(2)
    %h = histogram(y_hist, [-105:10:105]);
    %histogram(y_hist, [-105:10:105])
    %axis([-100 100 0 max(h)])

    %figure(3)
%     y_hist = y_hist';
%     plot(y_hist, (2*rand(length(y_hist),1)-1).*exp(-(y_hist/200).^2)/exp(1), '.')
%     axis([-10 10 -0.4 0.4])

%tiledlayout(2,1)

%nexttile
plot(z_hist,y_hist,'.')
xlabel('y')
ylabel('z')
%axis([-200 200 -(5*pi - pi/2) (5*pi - pi/2)])

% nexttile
% histogram(y_hist, 'normalization','pdf')
% hold on
% [kde, xi] = ksdensity(y_hist);
% plot(xi,kde,'linewidth',2)
% hold off



toc