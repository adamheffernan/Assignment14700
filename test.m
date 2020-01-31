clear all 
close all 

m_0 = 9.10938e-31;
m = 0.26*m_0;
T = 300;
k = 1.38064e-23;
v_th = sqrt((2*k*T)/m); 
mean = 0.2e-12;
x = v_th*mean; 

width = 100e-9;
len = 200e-9;
par = 1000;
pop = 50; 
iter = 500;
pos_velo = zeros(par,4);
temp = zeros(iter,1);
traj = zeros(iter, pop*2);
step = width/v_th/100;
for j = 1:par
    ang = rand*2*pi;
    pos_velo(j,:) = [len*rand width*rand v_th*cos(ang) v_th*sin(ang)];
end

for j = 1:iter
    
    pos_velo(:,1:2) = pos_velo(:,1:2) + step*pos_velo(:,3:4);
    
    i = pos_velo(:,1) > len;
    pos_velo(i,1) = pos_velo(i,1) - len;
    
    i = pos_velo(:,1) < 0;
    pos_velo(i,1) = pos_velo(i,1) + len;
    
    i = pos_velo(:,2) > width;
    pos_velo(i,2) = 2*width - pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
    
    i = pos_velo(:,2) < 0;
    pos_velo(i,2) = -pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
   x_vec(1:iter,j) = pos_velo(1:iter,1);
   y_vec(1:iter,j)= pos_velo(1:iter,2);
    
    plot(x_vec(1:iter,:)./1e-9, y_vec(1:iter,:)./1e-9, '.');
 
    
    figure(1);
    
        hold on;
        pause(0.05);
end
       
scat = 1 - exp(-step/mean);
v_boltz = makedist('Normal','mu',0,'sigma',sqrt(k*T/m));

for j = 1:par
    ang = rand*2*pi;
    pos_velo(j,:) = [len*rand width*rand random(v_boltz) random(v_boltz)];
     for i=1:pop
        traj(j, (i):(i+1)) = pos_velo(i, 1:2);  
     end 
end
for j = 1:iter
    pos_velo(:,1:2) = pos_velo(:,1:2) + step*pos_velo(:,3:4);
    
    i = pos_velo(:,1) > len;
    pos_velo(i,1) = pos_velo(i,1) - len;
    
    i = pos_velo(:,1) < 0;
    pos_velo(i,1) = pos_velo(i,1) + len;
    
    i = pos_velo(:,2) > width;
    pos_velo(i,2) = 2*width - pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
    
    i = pos_velo(:,2) < 0;
    pos_velo(i,2) = -pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
    
    i = rand(pop,1) < scat; 
    pos_velo(i,3:4) = random(v_boltz,[sum(i),2]);
  
    
        figure(3);
        subplot(2,1,1);
        plot(pos_velo(1:pop,1), pos_velo(1:pop,2), '.');
        axis([0 len 0 width]);
        title (['Trajectories for ', num2str(pop), ' of ', num2str(par), ' Electrons with Fixed Velocity'])
        xlabel 'x position (m)'
        ylabel 'y position (m)'
        subplot(2,1,2);
        plot(step*(0:j-1), temp(1:j))
        axis([0 step*iter 0 600])
        title 'Semiconductor temp'
        xlabel 'Time (s)' 
        ylabel 'temp (K)' 
        pause(0.02);
        figure(4);
        hold on;
        plot(traj(:,j)./1e-9, traj(:,j+1)./1e-9, '.');
end
%figure(4);
%hold on;
%for j=1:pop
 %plot(traj(:,j)./1e-9, traj(:,j+1)./1e-9, '.');
%end