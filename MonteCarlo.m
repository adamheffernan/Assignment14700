%Adam Heffernan 100977570 Assignment 1 4700. Completed on 2/1/2020. 

clear all 
close all 

m_0 = 9.10938e-31; %Mass of Electron
m = 0.26*m_0;%Effective mass in silicon
T = 300;%temperature in Kelvin 
k = 1.38064e-23;%Boltzmans Constant
v_th = sqrt((2*k*T)/m); %Thermal velocity of particles
mean = 0.2e-12;%Mean value 
width = 100e-9;%Width of Wafer
len = 200e-9;%Length of the channel 
par = 10000;%Particles
pop = 5; %Mapping population
iter = 1000;%Iterations
pos_velo = zeros(par,4);%Initializing position velocity Matrix
temp = zeros(iter,1);%Initializing position velocity Matrix
traj = zeros(iter, 2*pop);%Initializing trajectory Matrix
step = width/v_th/100;%Calculating step size for the movement of particles
for j = 1:par
        ang = rand*2*pi;
        pos_velo(j,:) = [len*rand width*rand v_th*cos(ang) v_th*sin(ang)];
end

for j = 1:iter
        %
        %Updating the position velocity vector at each iteration. Using step
        %size * position velocity vector of each particle at each iteration in time plus the
        %previous X and Y locations of the particle respectively. 
        pos_velo(:,1:2) = pos_velo(:,1:2) + step*pos_velo(:,3:4);
        %
        %Left Bound Check
        i = pos_velo(:,1) > len;
        pos_velo(i,1) = pos_velo(i,1) - len;
        %
        %Right Bound Check
        %
        i = pos_velo(:,1) < 0;
        pos_velo(i,1) = pos_velo(i,1) + len;
        %
        %Top Bound Check
        i = pos_velo(:,2) > width;
        pos_velo(i,2) = 2*width - pos_velo(i,2);
        pos_velo(i,4) = -pos_velo(i,4);
        %
        %Bottom Bound Check
        i = pos_velo(:,2) < 0;
        pos_velo(i,2) = -pos_velo(i,2);
        pos_velo(i,4) = -pos_velo(i,4);
        %
        %Storing the Temperature value of the semiconductor at each iteration
        temp(j) = (sum(pos_velo(:,3).^2) + sum(pos_velo(:,4).^2))*m/k/2/par;
        %Storing the Trajectory value of the semiconductor at each iteration
        traj(1:iter,2*j:2*j+1) = [pos_velo(1:iter,1) pos_velo(1:iter,2)];
        figure(1);
        subplot(2,1,1)
        plot(pos_velo(1:pop,1)./1e-9, pos_velo(1:pop,2)./1e-9, 'o');
        axis([0 len/1e-9 0 width/1e-9]);
        title('Simulation of Electrons in Silicon Crystal')
        xlabel('(nm)')
        ylabel('(nm)')
        subplot(2,1,2)
        plot(step*(0:j-1), temp(1:j),'Color',[0.1 0.1 0.1]);
        axis([0 step*iter 200 400]);
        title('Temperature of Electrons in Silicon Crystal')
        xlabel('Time(s)')
        ylabel('Temperature(K)')
        hold on;
end
    

for i = 1:pop
        color = [rand rand rand];
        figure(2)  
        plot(traj(i,2:2:end)./1e-9, traj(i,1:2:end-1)./1e-9, '.', 'Color',color);
        axis([0 len/1e-9 0 width/1e-9]);
        title('Trajectories of Electrons in Silicon Crystal')
        xlabel('(nm)')
        ylabel('(nm)')
        hold on;   
        
end
    scat = 1 - exp(-step/mean);
    v_boltz = makedist('Normal','mu',0,'sigma',sqrt(k*T/m));

for j = 1:par
        ang = rand*2*pi;
        pos_velo(j,:) = [len*rand width*rand random(v_boltz) random(v_boltz)];
end
   
    
for j = 1:iter
    
    pos_velo(:,1:2) = pos_velo(:,1:2) + step*pos_velo(:,3:4);
    
    traj(1:iter,2*j:2*j+1) = [pos_velo(1:iter,1) pos_velo(1:iter,2)];
    t = pos_velo(:,1) > len;
    pos_velo(t,1) = pos_velo(t,1) - len;
    
    i = pos_velo(:,1) < 0;
    pos_velo(i,1) = pos_velo(i,1) + len;
    
    i = pos_velo(:,2) > width;
    pos_velo(i,2) = 2*width - pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
    
    i = pos_velo(:,2) < 0;
    pos_velo(i,2) = -pos_velo(i,2);
    pos_velo(i,4) = -pos_velo(i,4);
    
    i= rand(par, 1) < scat;
    pos_velo(i,3:4) = random(v_boltz, [sum(i),2]);
    %Storing the Temperature value of the semiconductor at each iteration
    temp(j) = (sum(pos_velo(:,3).^2) + sum(pos_velo(:,4).^2))*m/k/2/par;
    %Storing the Trajectory value of the semiconductor at each iteration
    traj(1:iter,2*j:2*j+1) = [pos_velo(1:iter,1) pos_velo(1:iter,2)]; 
    figure(3)
    subplot(3,1,1)
    plot(pos_velo(1:pop,1)./1e-9, pos_velo(1:pop,2)./1e-9, 'o')
    axis([0 len/1e-9 0 width/1e-9])
    title('Simulation of Electrons in Silicon Crystal')
    xlabel('(nm)')
    ylabel('(nm)')
    subplot(3,1,2)
    plot(step*(0:j-1), temp(1:j),'Color',[0.1 0.1 0.1])
    axis([0 step*iter 200 400]);
    title('Temperature of Electrons in Silicon Crystal')
    xlabel('Time(s)')
    ylabel('Temperature(K)')
    hold on;
end
%Outputting Average Velocity
    velocity_average= sqrt(sum((pos_velo(:,3).^2))/par + sum(pos_velo(:,4).^2)/par);
    sprintf('Average Electron Velocities in Silicon Crystal are %s km/s',num2str(velocity_average/10e3));
%Output the average temperature 
    av_temp = sum(temp)/iter;
    sprintf('Average Temperature in the Crystal is %s K and the Average Electron Velocities in Silicon Crystal are %s m/s',num2str(av_temp),num2str(velocity_average))
for i = 1:pop
        
        
        color = [rand rand rand];
        figure(4);
        plot(traj(i,2:2:end)./1e-9, traj(i,1:2:end-1)./1e-9, '.', 'Color',color)
        axis([0 len/1e-9 0 width/1e-9])
        title('Trajectories of Electrons in Silicon Crystal')
        xlabel('(nm)')
        ylabel('(nm)')
        hold on
        
end
    figure(3);
    subplot(3,1,3);
    velocity = sqrt(pos_velo(:,3).^2 + pos_velo(:,4).^2);
    histogram(velocity)
    title('Histogram of Electron Velocities in Silicon Crystal')
    xlabel('Velocity (m/s)')
    ylabel('Frequency')
    boxes = 1e-9.*[80 120 0 40; 80 120 60 100];
    boxes_specular = [0 1];
for j = 1:par
        ang = rand*2*pi;
        pos_velo(j,:) = [len*rand width*rand random(v_boltz) random(v_boltz)];
            while(in_box(pos_velo(j,1:2), boxes))
                    pos_velo(j,1:2) = [len*rand width*rand];
            end
end
