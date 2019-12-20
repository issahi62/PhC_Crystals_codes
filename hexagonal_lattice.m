%% INITIALIZE 
clear 
close all 
clc 

%% primitive lattice 
a = 1; 

a1 = a.*[1 0 0 ]; 
a2 = [a/2, sqrt(3)/2*a, 0];
a3 = [0, 0 a]; 

%% basis
b1= [0 0 0]; 

%% Figure
figure
hold on 

%%Draw single vector
quiver3(b1(1), b1(2), b1(3), a1(1), a1(2), a1(3), 'g', 'Linewidth', 2); 
quiver3(b1(1), b1(2), b1(3), a2(1), a2(2), a2(3), 'r', 'Linewidth', 2);
quiver3(a2(1), a2(2), a2(3), a3(1), a3(2), a3(3), 'b', 'Linewidth', 2);
grid on

N = 0; 
for n1 = 0:2
    for n2 = 0:2 
        for n3 = 0:2
            N = N+1;
            R(N, :) = n1.*a1 + n2.*a2 + n3.*a3 +b1 ; 
        end 
    end 
end

scatter3(R(:, 1), R(:,2), R(:,3), 500, 'MarkerFaceColor', [0, .75, .75]); 