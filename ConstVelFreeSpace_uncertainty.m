tic

sample_size = 100000;
speed = 0.1;

set(0,'DefaultAxesFontSize',20)

    x = randn(1,sample_size);
    x = x/(max(x) - min(x));
    min_x = min(x);
    v=0*x;

    omega = 95;
    nu = 1/2.1;
    K = omega/2/pi;
    C = 1/3;

    a = rand(1,sample_size);
    b = rand(1,sample_size);

%     writerObj = VideoWriter('FreeSpace.avi');
%     writerObj.FrameRate = 10; 
%     open(writerObj);

iter = 10000;
pdx = zeros(1, iter);
pdv = pdx;
    
for i = 1:iter

    b = 4*b.*(1-b);
    a = mod(a + b*sqrt(2),1);

    v = C*(v + K*cos(2*pi*a).*sin(x).*exp(-nu*abs(v)));
    v_actual = v/50 + speed;
    x = x + v_actual;
 
    pdx_fit = fitdist(x', 'Normal');
    pdx(i) = pdx_fit.sigma;
    pdv(i) = fitdist(v_actual', 'Normal').sigma;

end

pdx_fft = 1./pdx;
plot(pdx, pdv,'.', pdx, pdx_fft,'.')

toc