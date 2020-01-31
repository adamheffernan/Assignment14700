clear all 
close all 

%%%%%%%%%% constants %%%%%%%%%%
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
pop = 5; 
iter = 100;

%%%%%%%%%% Struct %%%%%%%%%%
point.pos_x = 0;
point.pos_y = 0;
point.vel_x = 0;
point.vel_y = 0;
lines = repmat(point, pop, iter);

step = width/v_th/100;
scat = 1 - exp(-step/mean);
v_boltz = makedist('Normal','mu',0,'sigma',sqrt(k*T/m));
for i = 1:par
    ang = rand*2*pi;
    point.pos_x = len*rand;
    point.pos_y = width*rand;
    point.vel_x = v_th*cos(ang);
    point.vel_y = v_th*sin(ang);
    lines(i, 1) = point;
end

figure(1);
subplot(2,1,1);
axis([0 len 0 width]);
title (['Trajectories for ', num2str(pop), ' of ', num2str(par), ' Electrons with Fixed Velocity'])
xlabel 'x position (m)'
ylabel 'y position (m)'
hold on;

for i = 1:pop
    color = [rand, rand, rand];
    for j = 2:iter
        % Move particle in direction of its velocity
        lines(i,j).pos_x = lines(i,j-1).pos_x + step*lines(i,j-1).vel_x;
        lines(i,j).pos_y = lines(i,j-1).pos_y + step*lines(i,j-1).vel_y;
        lines(i,j).vel_x = lines(i,j-1).vel_x;
        lines(i,j).vel_y = lines(i,j-1).vel_y;
        
        % Check if particle touches top
        if lines(i,j).pos_x > len
            lines(i,j).pos_x = lines(i,j).pos_x - len;
        % Check if particle touches bottom
        elseif lines(i,j).pos_x < 0
            lines(i,j).pos_x = lines(i,j).pos_x + len;
        end
        
        % Check if particle touches right
        if lines(i,j).pos_y > width
            lines(i,j).pos_y = 2*width - lines(i,j).pos_y;
            lines(i,j).vel_y = -lines(i,j).vel_y;
            
        % Check if particle touches left
        elseif lines(i,j).pos_y < 0
            lines(i,j).pos_y = -lines(i,j).pos_y;
            lines(i,j).vel_y = -lines(i,j).vel_y;
            
        
        end
        
        %[plotting
        plot(lines(i,j).pos_x, lines(i,j).pos_y, 'Marker', '.', 'Color', color);
        pause(0.00001);
    end
end



for i = 1:par
    ang = rand*2*pi;
    point.pos_x = len*rand;
    point.pos_y = width*rand;
    point.vel_x = random(v_boltz);
    point.vel_y = random(v_boltz);
    lines(i, 1) = point;
end

figure(2);
subplot(2,1,1);
axis([0 len 0 width]);
title (['Trajectories for ', num2str(pop), ' of ', num2str(par), ' Electrons with Fixed Velocity'])
xlabel 'x position (m)'
ylabel 'y position (m)'
hold on;

for i = 1:pop
    color = [rand, rand, rand];
    for j = 2:iter
        % Move particle in direction of its velocity
        lines(i,j).pos_x = lines(i,j-1).pos_x + step*lines(i,j-1).vel_x;
        lines(i,j).pos_y = lines(i,j-1).pos_y + step*lines(i,j-1).vel_y;
        lines(i,j).vel_x = lines(i,j-1).vel_x;
        lines(i,j).vel_y = lines(i,j-1).vel_y;
        t= rand(pop,1) < scat ;
        lines(t,t).vel_x = random(v_boltz,sum(t));
        lines(t,t).vel_y = random(v_boltz,sum(t));
        
        % Check if particle touches top
        if lines(i,j).pos_x > len
            lines(i,j).pos_x = lines(i,j).pos_x - len;
        % Check if particle touches bottom
        elseif lines(i,j).pos_x < 0
            lines(i,j).pos_x = lines(i,j).pos_x + len;
        end
        
        % Check if particle touches right
        if lines(i,j).pos_y > width
            lines(i,j).pos_y = 2*width - lines(i,j).pos_y;
            lines(i,j).vel_y = -lines(i,j).vel_y;
            
        % Check if particle touches left
        elseif lines(i,j).pos_y < 0
            lines(i,j).pos_y = -lines(i,j).pos_y;
            lines(i,j).vel_y = -lines(i,j).vel_y;
            
        
        end
        
        %[plotting
        plot(lines(i,j).pos_x, lines(i,j).pos_y, 'Marker', '.', 'Color', color);
        pause(0.00000000001);
    end
end



