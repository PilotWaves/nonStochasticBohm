tic

sample_size = 100001;
y_hist = 0;
z_hist = 0;
count = 1;

for i = 1:sample_size

    N = 1001;
    N_half = (N-1)/2 + 1;

    x = 0:N-1;
    vy = 0*x;
    y = vy;
    a = vy;
    b = vy;

    vz = vy;

    omega = 95;
    nu = 1/2.1;
    K = omega/2/pi;
    C = 1/3;

    a(1) = 0;
    b(1) = 0.1;

    y(1) = randn;
    z = y;

for n = 1:N_half-1
    
    b(n+1) = 4*b(n)*(1-b(n));
    a(n+1) = mod(a(n) + b(n+1)*sqrt(2),1);

    vy(n+1) = C*(vy(n) + K*cos(2*pi*a(n+1))*sin(2*pi*y(n))*exp(-nu*abs(vy(n))));
    y(n+1) = y(n) + vy(n+1);
    vz(n+1) = C*(vz(n) + K*cos(2*pi*a(n+1))*cos(2*pi*z(n))*exp(-nu*abs(vy(n))));
    z(n+1) = z(n) + vz(n+1);

end

    if y(N_half) > -5 && y(N_half) < 5

        
        for n = N_half:N-1

            b(n+1) = 4*b(n)*(1-b(n));
            a(n+1) = mod(a(n) + b(n+1)*sqrt(2),1);

            vy(n+1) = C*(vy(n) + K*cos(2*pi*a(n+1))*sin(y(n))*exp(-nu*abs(vy(n))));
            y(n+1) = y(n) + vy(n+1);
            vz(n+1) = C*(vz(n) + K*cos(2*pi*a(n+1))*cos(2*pi*z(n))*exp(-nu*abs(vy(n))));
            z(n+1) = z(n) + vz(n+1);
        end

    else

        x(N_half + 1: N) = fliplr(x(1:N_half-1));
        y(N_half + 1) = 2*y(N_half) - y(N_half-1);
        b(N_half + 1) = b(N_half);
        a(N_half + 1) = a(N_half);
        vy(N_half + 1) = vy(N_half - 1);
        
        for n = N_half:N-1

            b(n+1) = 4*b(n)*(1-b(n));
            a(n+1) = mod(a(n) + b(n+1)*sqrt(2),1);

            vy(n+1) = C*(vy(n) + K*cos(2*pi*a(n+1))*sin(2*pi*(y(n)-y(N_half)))*exp(-nu*abs(vy(n))));
            y(n+1) = y(n) + vy(n+1);
            vz(n+1) = C*(vz(n) + K*cos(2*pi*a(n+1))*cos(2*pi*(z(n)-z(N_half)))*exp(-nu*abs(vy(n))));
            z(n+1) = z(n) + vz(n+1);
        end

    end

    %plot(y, x)
    %hold on

    if y(N_half) > -5 && y(N_half) < 5
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

figure(1)
plot(y_hist,z_hist,'.')

toc