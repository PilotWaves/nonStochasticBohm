tic

slits = 3;

sample_size = 100001;
y_hist = zeros(sample_size,1);
count = 1;

for i = 1:sample_size

    N = 501;

    x = 0:N-1;
    v = 0*x;
    y = v;
    a = v;
    b = v;

    nu = 1/2.1;
    K = 15;
    C = 1/3;

    a(1) = 0;
    b(1) = 0.1;

    y(1) = 6e-3*i/sample_size - 3e-3;

for n = 1:N-1
    
    b(n+1) = 4*b(n)*(1-b(n));
    a(n+1) = mod(a(n) + b(n+1)*sqrt(2),1);

    v(n+1) = C*(v(n) + K*cos(2*pi*a(n+1))*sin(pi*y(n)/slits)*exp(-nu*abs(v(n))));
    y(n+1) = y(n) + v(n+1);

end

%         plot(x,y, 'linewidth',0.5)
%         hold on
    y_hist(i) = y(end);

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
    plot(y_hist, rand(length(y_hist),1), '.')
    axis([-10 10 0 1])
toc