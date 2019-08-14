discrete = 1;
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

H = kine + V*1;
U = expm(-1i*H*1e1);
init  = zeros(N,1);
init((N-1)/32*15) = 1;
state = init;
while 1
    
    
    %p_state = state;
    state = U*state;
    plot(fs/k_R*2*pi,sq(state))
    drawnow;
end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'conv.gif';
dt = 5/pi;
t=0;
init  = zeros(N,1);
init((N-1)/32*15) = 1;
state = init;
for n = 1:2000
    % Draw plot for y = x.^n
    t = t+dt;
    state = U*state;
    if mod(n,20)==0
        plot(fs/k_R*2*pi,sq(state))
        xlabel('k/k_R')
        title(strcat('t=',num2str(t),'ms'))
        drawnow;
          % Capture the plot as an image 
          frame = getframe(h); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if n == 20 
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
              imwrite(imind,cm,filename,'gif','WriteMode','append'); 
          end 
    end
end



