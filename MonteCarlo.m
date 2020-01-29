clear all 
close all 

m_0 = 9.10938e-31;
m = 0.26*m_0;
T = 300;
k = 1.38064e-23;
v_th = sqrt((2*k*T)/m); 

x = v_th*0.2e-12; 

ht = 100e-9;
ln = 200e-9;
par = 1000;
pop = 10; 
iter = 1000;
state = zeros(par,4);
temp = zeros(iter,1);
step = ht/v_th/100;
for j = 1:par
    ang = rand*2*pi;
    state(j,:) = [ln*rand ht*rand v_th*cos(ang) v_th*sin(ang)];
end

for j = 1:iter
    state(:,1:2) = state(:,1:2) + step*state(:,3:4);
    i = state(:,1) > ln;
    state(i,1) = state(i,1) - ln;
    
    i = state(:,1) < 0;
    state(i,1) = state(i,1) + ln;
    
    i = state(:,2) > ht;
    state(i,2) = 2*ht - state(i,2);
    state(i,4) = -state(i,4);
    
    i = state(:,2) < 0;
    state(i,2) = -state(i,2);
    state(i,4) = -state(i,4);
    
     figure(1);
        subplot(2,1,1);
        
        plot(state(1:pop,1)./1e-9, state(1:pop,2)./1e-9, 'o');
        axis([0 ln/1e-9 0 ht/1e-9]);
        title(sprintf('Trajectories for %d of %d Electrons with Fixed Velocity (Part 1)',...
        pop,par));
        xlabel 'x (nm)'
        ylabel 'y (nm)'
        
            subplot(2,1,2);
            
            plot(step*(0:j-1), temp(1:j))
            axis([0 step*iter 0 600])
            title 'Semiconductor temp'
            xlabel 'Time (s)' 
            ylabel 'temp (K)' 
            pause(1);
end