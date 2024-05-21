%%

tic

sample_size = 100001;
speed = 0.1;

%sample_size = 1000;
%speed = 0;

wall = 10;
wall_height = speed;
wall_width = speed;
boundary = 5;

N = 500;



x = zeros(4,N,sample_size);

v = x;
a = x;
b = x;
speed = ones(4,N,sample_size)*speed;

C = 1/3;
nu = 1/2.1;
K_vec = [0, 10, 100,1000];

for k = 4:-1:1

x(k,1,:) = randn(1,sample_size);
%x(k,1,:) = -boundary + boundary*2*rand(1,sample_size);
a(k,1,:) = rand(1,sample_size);
b(k,1,:) = rand(1,sample_size);

x(k,1,:) = x(k,1,:)/(max(x(k,1,:)) - min(x(k,1,:)));

K = K_vec(k);

for i = 1:N-1
   for j = 1:sample_size

    b(k,i+1,j) = 4*b(k,i,j)*(1-b(k,i,j));
    a(k,i+1,j) = mod(a(k,i,j) + b(k,i+1,j)*sqrt(2),1);

    v(k,i+1,j) = C*(v(k,i,j) + K*cos(2*pi*a(k,i+1,j))*sin(x(k,i,j))*exp(-nu*abs(v(k,i,j))));
    x(k,i+1,j) = x(k,i,j) + v(k,i+1,j)/50 + speed(k,i,j);
    speed(k,i+1,j) = speed(k,i,j);

    if abs(x(k,i+1,j)) > boundary
            v(k,i+1,j) = -v(k,i+1,j);
            x(k,i+1,j) = sign(x(k,i+1,j))*2*boundary - x(k,i+1,j);
            speed(k,i+1,j) = -speed(k,i+1,j);
    end

   end

   %[k, i]

end


end

save("ParticleInABox_Vars.mat", "x")



toc

%%

tic
set(0,'DefaultAxesFontSize',20)

writerObj = VideoWriter('ParticleInABox_test_vel.avi');
writerObj.FrameRate = 10; 
open(writerObj);

%box_axis = [-boundary boundary 0 7; -boundary boundary 0 3; -boundary boundary 0 0.75; -boundary boundary 0 0.25];

sample_size = length(x(1,1,:));
N = length(x(1,:,1));

%x_kde = zeros(sample_size,1);
v_kde = zeros(sample_size,1);
max_kde = [7,7,7,7];

for i = 1:N
    tiledlayout(4,1,'TileSpacing','compact')

    for k = 1:4
    nexttile
    %histogram(x(k,i,:), 'Normalization','pdf','edgecolor','none')

    x_kde(:) = x(k,i,:);
    %v_kde(:) = v(k,i,:);

    [kde, xi] = ksdensity(x_kde);
    plot(xi,kde,'linewidth',2)
    
    if max(kde) < max_kde(k)/5 && k == 2
        max_kde(k) = max(kde)*1.75;
    elseif max(kde) < max_kde(k)/5 && k == 3
        max_kde(k) = max(kde)*1.75;
    elseif max(kde) < max_kde(k)/5 && k == 4
        max_kde(k) = max(kde);
    end

    axis([-5 5 0 max_kde(k)])
    xticks([])

    % 
    % nexttile

    % [kde, xi] = ksdensity(v_kde);
    % plot(xi,kde,'linewidth',2)
    %axis(box_axis(k,:))

    end

    frame = getframe(gcf);
    writeVideo(writerObj,frame);

end

close(writerObj)
toc

%%

t_layout = tiledlayout(2,2,'TileSpacing','Compact');
for k = 1:4
nexttile
for i = 350:25:500
x_kde(:) = x(k,i,:); [kde, xi] = ksdensity(x_kde);
plot3(xi, i*ones(1,length(xi)), kde, 'linewidth', 2)
hold on
end
xticks([])
yticks([])
hold off
end
xlabel(t_layout, 'x','fontsize',20)
ylabel(t_layout, 'PDF','fontsize',20)
