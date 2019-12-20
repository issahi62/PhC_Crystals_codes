% **************
%% INITIALIZE 
%****************
clc 
close all 
clear 
% **************
%% DASHBOARD 
%****************
n = 1; 
er1 = 1; ur1=1; er2=1; ur2 =1; 
theta =pi/2; 
phi = pi; 
pte = 1; 
ptm =1;
lambd0=633; 
K0 = 2*pi/lambd0;

% **********************************
%% Calculate Transverse Wave Vectors
%***********************************
Kx = n.*sin(theta).*cos(phi); 
Ky = n.*sin(theta).*sin(phi); 

% **********************************
%% Calculate Free Space parameters
%***********************************
Kz = sqrt(1-Kx.^2-Ky.^2); 
Q = [Kx.*Ky, 1-Kx.^2; 
    Ky.^2-1 -Kx.*Ky]; 
I =eye(size(Q,2));
Gamma = 1j.*Kz.*I;
V = Q.*inv(Gamma); 
S11_const = ones(size(Q, 2));
S21_const = ones(size(Q, 2));

% **********************************
%% INITIALIZE GLOBAL SCATTERING MATRIX
%***********************************
SG11 = 0.*I;
SG22 = 0.*I;
SG12 =I;
SG21 =I; 
SG = [SG11, SG12;
            SG21, SG22]; 
% **************
%% LAYERS
%****************
ER  = [2.5,  3.5 ,2]; %layers  
UR =  [1, 1, 1]; %permeability 
L =   [0.2 0.3 0.5]; % length of propagation 

% **********************************
%% INITIALIZE Loop
%***********************************
for i = 1:length(ER) 
    
% **********************************
%% Parameters for layers
%***********************************
    Kz(i) = sqrt(UR(i)*ER(i)-Kx.^2 -Ky.*2); 
    Q(:, :, i) = 1/UR(i).*[Kx.*Ky, UR(i)*ER(i)-Kx.^2; 
    Ky.^2-UR(i)*ER(i) -Kx.*Ky];
    Gamma(:, :, i) =1j.*Kz(i).*I;
    V(:, :, i) = Q(:, :, i).*inv(Gamma(:, :, i));
    
%**********************************
%% Scattering Matrix for Layer i 
%***********************************

A(:, :, i) = I + (V(:, :, i)^-1)*V(:, :, i);
B(:, :, i) = I - (V(:, :, i)^-1)*V(:, :, i);

X(i) = exp(-lambd0*K0*L(i)); 

S11(:, :, i) = (A(:, :, i)-X(i)*B(:,:,i).*inv(A(:, :, i)).*X(i).*...
    B(:, :, i))^-1 .* (X(i)*B(:,:,i).*inv(A(:, :, i)).*X(i).*...
    A(:, :, i) - B(:, :, i));
S21(:, :, i) = (A(:, :, i)-X(i)*B(:,:,i).*inv(A(:, :, i)).*X(i).*...
    B(:, :, i))^-1 .* X(i)*(A(:, :, i)- B(:, :, i).*inv(A(:, :, i)).*...
    B(:, :, i)); 

S11AB(:, :, i) = SG11 + SG21.*inv(I - S11(:, :, i).*SG22).*(S11(:, :, i).*SG21);
S21AB(:, :, i) = SG21.*inv(I - S11(:, :, i).*SG22).*SG21;
%S11_const = S11_const.*S11(:, :,i);
%S21_const = S21_const.*S21(:, :, i); 

%**********************************
%% UPDATE SCATTERING
%***********************************
end 


% **************
%% TMM 
%****************

% G = [0 0 0 ur; 
%      0 0 -ur 0; 
%      0 eps 0 0; 
%      -eps 0 0 0];
%  
%  K = eye(size(G)); 
%  
%  [V, D] = eig(G);
%  D = sort(diag(D), 'descend'); 
%  V = sort(V); 
%  D = K.*D; 
%  
%V = [V(:, 1), V(:, 3), V(:, 2), V(:, 4)]; 
 
 