%Daniel Adea
%204999515
%Euler-Bernoulli Beam Bending
%This problem models a simply-supported aluminum beam that is subjected to
%a single point force

%clear workspace
clear all; clc; close all;

%set variables as given 
nodes = 20;
L = 1;
nodeDist = L/(nodes-1);
R = .013;
r = .011;
P = -2000;
E = 70e9;
I = pi/4*(R^4 - r^4);
d = 0.75;

%preallocate arrays for later use
points = linspace(0, L, nodes);
b = zeros(nodes, 1);
A = zeros(nodes, nodes);

%set array A

    %set the endpoints
A(1,1) = 1;
A(nodes, nodes) = 1;

    %use a for loop to set inner values following the pattern
for i = 2:nodes-1
    A(i,i) = -2;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
end

%set matrix b
    %set the boundary conditions
b(1) = 0;
b(nodes)= 0;

    %calculate the right side of the equation
for j = 2:nodes-1 
    pointDist = (j-1)*nodeDist;
    %use the correct equation depending on location in comparison to force
    if pointDist <= d
        b(j) = nodeDist^2*(-P*(L-d)*pointDist/L)/(I*E);
    else
        b(j) = nodeDist^2*(-P*d*(L-pointDist)/L)/(I*E);
    end
end

%calculate displacement
y = A\b;

%plot displacement
plot(points, y, 'o-');

%calculate theoretical maximum displacement
c = min(d, L-d);
ymax = P*c*(L^2 - c^2)^1.5/(9*sqrt(3)*E*I*L);




    
    