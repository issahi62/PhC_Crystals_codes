%% INITIALIZE

clc 
close all;
clear; 

%% Translation vectors
a1 = [ 1 0 0]; 
a2 = [ 0 1 0]; 
a3 = [0  0 1]; 

%% BAsis atoms
b1 = [0 0 0]; 
b2 = [0 0.5 0.5]; 
b3 = [0.5 0 0.5]; 
b4 = [0.5 0.5 0]; 

% 
figure 
hold on 
%% arrows 

quiver3(b1(1), b1(2), b1(3), a1(1), a1(2), a1(3), 'g', 'Linewidth', 2); 
quiver3(b1(1), b1(2), b1(3), a2(1), a2(2), a2(3), 'r', 'Linewidth', 2);
quiver3(b1(1), b1(2), b1(3), a3(1), a3(2), a3(3), 'b', 'Linewidth', 2);
grid on 
N=0; 
%% SUPERCELL
for n1 = 0:2
    for n2 = 0:2 
        for n3 = 0:2 
            N= N+1;
            R(N, :) = (n1.*a1+ n2.*a2 + n3.*a3)+b1;
            
            N= N+1;
            R(N, :) = n1.*a1+n2.*a2+n3.*a3 +b2; 
            
            N= N+1;
            R(N, :) = n1.*a1+n2.*a2+n3.*a3 +b3;
            
            N= N+1;
            R(N, :) = n1.*a1+n2.*a2+n3.*a3 +b4;
        end 
    end 
end 

%% Plot supercell
scatter3(R(:, 1), R(:,2), R(:, 3), 500, 'MarkerFaceColor', [0, .56, .56]); 

            