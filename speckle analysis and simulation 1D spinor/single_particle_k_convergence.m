Rb_Mass = 1.42*10^-25;% mass of Rb87 in Kg
hbar = 1.05*10^-34;
raman_wavelength = 790*10^-9;
N = 2^10 + 1;
k_R = 2*pi/raman_wavelength;
dk = 2*k_R/(N-1);
dx = 1/N/dk;
X = ((0:1:N-1)-(N-1)/2).*dx;
K = ((0:1:N-1)-(N-1)/2).*dk;

L = dx*(N-1);


kine = diag((K./k_R).^2*3.7);

harmonic = 1/2*Rb_Mass*(2*pi*100)^2.*X.^2/hbar/2/pi/1000;
harmonic = harmonic - mean(harmonic);
harmonicK = fourier_transform(harmonic,N,dx)/L;
rev = fliplr(harmonicK);

V = zeros(N,N);

mid = (N-1)/2 + 1;

V(:,mid) = rev;


for i=1:mid-1
    V(1:N-i,mid-i) = rev(1+i:end);
    V(1+i:end,mid+i) = rev(1:N-i);
end


dt = 0.05;
H = kine + V;
U = expm(-1i*H*dt);
% init  = zeros(N,1);
% init((N-1)/32*15) = 1;
% state = init;
% while 1
%     
%     
%     %p_state = state;
%     state = U*state;
%     plot(fs/k_R*2*pi,sq(state))
%     drawnow;
% end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'convH100Hz.gif';
dtms = dt/2/pi;
t=0;
init  = zeros(N,1);
init((N-1)/8*3) = 1;
state = init;
for n = 1:5000
    % Draw plot for y = x.^n
    t = t+dtms;
    state = U*state;
    if mod(n,10)==0
        plot(K./k_R,sq(state))
        xlabel('k/k_R')
        title(strcat('t=',num2str(t),'ms'))
        drawnow;
          % Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if n == 2 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 2^(13-discrete)+1;
fs = ((0:1:N-1)-(N-1)/2).*2^discrete*deltaf;
P = zeros(N,Nx);
for i=1:N
    P(i,1+2^discrete*(i-1)) = 1;
end
%hbar^2k_R^2/2/m~2pi*3.7kHz, scale it to 3.7
%mean(v(x))=2pi*1kHz, scale it to 1
i=8;
filepath = strcat('speckle bench test data/numerical_speckle/13/inten_', num2str(i),'.mat');
speckle = load(filepath);
speckle = speckle.inten;

speckle = speckle/mean(mean(speckle));%average intensity about 1 in simulation units


speckle_k = mean(fourier_transform(speckle,Nx,deltax))/L;

harmonic = 1/2*Rb_Mass*(2*pi*100)^2.*(a_0*X).^2/hbar/2/pi/1000;
harmonic = harmonic - mean(harmonic);
harmonicK = fourier_transform(harmonic,Nx,deltax)/L;
rev = P*fliplr(harmonicK)';

rev = P*fliplr(speckle_k)';
%plot(f./k_R*2*pi,sq(speckle_k))
% rev =  zeros(N,1);
% rev((N-1)/4:(N-1)/2+1) = linspace(0,2,(N-1)/4+2);
% rev = (rev + flipud(rev))/2;
% rev((N-1)/2+1) = 0;
kine = diag((fs./k_R*2*pi).^2*3.7);

V = zeros(N,N);

mid = (N-1)/2 + 1;

V(:,mid) = rev;


for i=1:mid-1
    V(1:N-i,mid-i) = rev(1+i:end);
    V(1+i:end,mid+i) = rev(1:N-i);
end

H = kine + V*20;
%U = eye(N)-1i*H*1e-4;
U = expm(-1i*H*1e-1);
init  = zeros(N,1);
init(300) = 1;

difference = 1;
state = init;

% while difference>1e-10
%     p_state = state;
%     state = U*state;
%     difference = max(abs(sq(p_state) - sq(state)))
%     %sum(sq(state))
%     plot(fs/k_R*2*pi,sq(state))
%     drawnow;
% end

dt = 1e-1;
H = kine + V;
U = expm(-1i*H*dt);
% init  = zeros(N,1);
% init((N-1)/32*15) = 1;
% state = init;
% while 1
%     
%     
%     %p_state = state;
%     state = U*state;
%     plot(fs/k_R*2*pi,sq(state))
%     drawnow;
% end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'convH100Hz.gif';
dtms = dt/2/pi;
t=0;
init  = zeros(N,1);
init((N-1)/32*15) = 1;
state = init;
for n = 1:500
    % Draw plot for y = x.^n
    t = t+dtms;
    state = U*state;
    if mod(n,2)==0
        plot(fs/k_R*2*pi,sq(state))
        xlabel('k/k_R')
        title(strcat('t=',num2str(t),'ms'))
        drawnow;
          % Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if n == 2 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end



