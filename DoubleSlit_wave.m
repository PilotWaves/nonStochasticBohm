tic

sample_size = 10001;
y_hist = 0;
count = 1;
slit_pos = 30;
slit_width = 10;

tiledlayout(2,1)
nexttile

colormap(turbo)
plot([500 500], [-100 -30], 'k', [500 500], [-20 20], 'k', [500 500], [30 100], 'k', 'linewidth',5)
hold on

    N = 1001;
    N_half = 501;

    v = zeros(sample_size, N);
    y = v;
    x = v;
    a = v;
    b = v;

    omega = 95;
    nu = 1/2.1;
    K = omega/2/pi;
    C = 1/3;

    a(:,1) = 0;
    b(:,1) = 0.1;

    x(:,1) = 1;
    y(:,1) = 6e-3*[1:sample_size]/sample_size - 3e-3;
    


% [h, centers] = hist3([x(:,1), y(:,1)],'CDataMode','auto', 'Nbins', [50 50]);
% centers1 = centers{1,1};
% centers2 = centers{1,2};
% surf(centers1, centers2, sqrt(h'), 'FaceAlpha',0.4)
% view(2)
% xlim([0 N-1])
% ylim([-100 100])
% shading interp
% hold on

for n = 1:N_half-1
    
    x(:,n+1) = n+1;
    b(:,n+1) = 4*b(:,n).*(1-b(:,n));
    a(:,n+1) = mod(a(:,n) + b(:,n+1)*sqrt(2),1);

    v(:,n+1) = C*(v(:,n) + K*cos(2*pi*a(:,n+1)).*sin(y(:,n)).*exp(-nu*abs(v(:,n))));
    y(:,n+1) = y(:,n) + v(:,n+1);

% [h, centers] = hist3([x(:,n+1), y(:,n+1)],'CDataMode','auto', 'Nbins', [50 50]);
% centers1 = centers{1,1};
% centers2 = centers{1,2};
% surf(centers1, centers2, sqrt(h'), 'FaceAlpha',0.4)
% view(2)
% xlim([0 N-1])
% ylim([-100 100])
% shading interp
% hold on

end


for n = N_half:N-1
    X = 0;
    Y = 0;
    count = 0;
    for i = 1:sample_size

if y(i,N_half) > 60
    traj_bounce = i;
end

    if (y(i, N_half) > -slit_pos && y(i, N_half) < -slit_pos + slit_width) || (y(i, N_half) > slit_pos - slit_width && y(i, N_half) < slit_pos)


%         if y(N_half) > -slit_pos && y(N_half) < -slit_pos + slit_width
%             center = -slit_pos + slit_width/2;
%         end
% 
%         if y(N_half) > slit_pos - slit_width && y(N_half) < slit_pos
%             center = slit_pos - slit_width/2;
%         end


        traj_example = i;
            x(i,n+1) = n+1;
            b(i,n+1) = 4*b(i,n)*(1-b(i,n));
            a(i,n+1) = mod(a(i,n) + b(i,n+1)*sqrt(2),1);

            v(i,n+1) = C*(v(i,n) + K*cos(2*pi*a(i,n+1))*sin(y(i,n))*exp(-nu*abs(v(i,n))));
            y(i,n+1) = y(i,n) + v(i,n+1);

%             count = count + 1;
%             X(count) = x(i,n+1);
%             Y(count) = y(i,n+1);

    else

        x(i,n+1) = x(i, n) - 1;
        y(i, N_half + 1) = 2*y(i, N_half) - y(i, N_half-1);
        b(i, N_half + 1) = b(i, N_half);
        a(i, N_half + 1) = a(i, N_half);
        v(i, N_half + 1) = v(i, N_half - 1);

            b(i,n+1) = 4*b(i,n)*(1-b(i,n));
            a(i,n+1) = mod(a(i,n) + b(i,n+1)*sqrt(2),1);

            v(i,n+1) = C*(v(i,n) + K*cos(2*pi*a(i,n+1))*sin(y(i,n))*exp(-nu*abs(v(i,n))));
            y(i,n+1) = y(i,n) + v(i,n+1);
    end

    end

%[h, centers] = hist3([X', Y'],'CDataMode','auto', 'Nbins', [2 50]);
[h, centers] = hist3([x(:,n+1), y(:,n+1)],'CDataMode','auto', 'Nbins', [50 50]);
centers1 = centers{1,1};
centers2 = centers{1,2};
surf(centers1, centers2, sqrt(h'), 'FaceAlpha',0.4)
view(2)
xlim([0 N-1])
ylim([-100 100])
%pause(0.1)
shading interp
hold on

end
    hold off

    nexttile
    plot([500 500], [-100 -30], 'k', [500 500], [-20 20], 'k', [500 500], [30 100], 'k', 'linewidth',5)
hold on
plot(x(traj_bounce,:),y(traj_bounce,:),x(traj_example,:),y(traj_example,:),'linewidth',2)
hold off
xlim([0 N-1])
ylim([-100 100])

toc