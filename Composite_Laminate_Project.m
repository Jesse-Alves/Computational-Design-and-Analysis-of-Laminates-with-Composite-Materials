%% ============= Project of a Composite Laminate ==================
% Health Tech Master
% Date: 16/01/2023
% ----------------------------- Jesse Alves
clear all; clc;
%% ==================================================
%% ======== Calculate the % of folds in each direction ========
%% ==================================================
disp('============================================')
disp('===============  QUESTION 1==================')
disp('============================================')
% Parameters
Vf = 0.6;
e = 0.2;
E1 = 85000; E2 = 5600; G12 = 2100; v12 = 0.34; 
sigmaRL = 1410; sigmaRT = 28; tRLT = 45;
Nx = 1000; Ny = 1000; Txy = 200;

% Computing ex ey and exy
ex = Nx/sigmaRL;
ey = Ny/sigmaRL;
exy = Txy/tRLT;

disp('Percentage of each plies')
percentage_0 = (ex)/(ex + ey + exy)
percentage_90 = (ey)/(ex + ey + exy)

% percent 45° = -45°
percentage_45 = (exy)/(2*(ex + ey + exy))
percentage_minus45 = percentage_45

Total_Percentage = percentage_0 + percentage_45 + percentage_minus45 + percentage_90

message1 = strcat('The percentage of 0° and 90° is', {' '}, num2str(100*percentage_0),'%');
msgbox(message1)

message1 = strcat('The percentage of 45° and -45° is', {' '}, num2str(100*percentage_45),'%');
msgbox(message1)

%% ================================================================
%% == Calculate the minimum thickness of laminate required to prevent breakage.
%% ================================================================
disp('============================================')
disp('===============  QUESTION 2==================')
disp('============================================')

%% Obtaining the number of plies

disp(['After adding the number of plies to get the percentage, the best ' ...
    'approximation was:'])
disp(' ')
disp(' ')

disp('The desired percentage:')
percentage_90 = percentage_90
percentage_45 = percentage_45

disp('The percentage obtained with:  theta = [90 -45 45 -45 45 -45 45  0 0 -45 45 -45 45 -45 45  90] ')
theta = [90 -45 -45 -45 45 45 45  0 0 45 45 45 -45 -45 -45  90];

p_90 = 2/length(theta)
p_45 = 6/length(theta)

% As the percentage of 0 and 90 are the same and the 45 and -45 are the
% same, just the 90 and 45 was calculated.


%% Matrix Q
v21 = (E2/E1)*v12;

Q11 = (E1)/(1 - v12*v21);
Q12 = (v21*E1)/(1-v12*v21); %Q12 = (v12*E2)/(1-v12*v21);
Q22 = (E2)/(1-v12*v21);
Q33 = 2*G12;

% The Q matrix
Q = [Q11 Q12 0;
          Q12 Q22 0;
           0     0   Q33];

%% Matrix Q bar

theta = [90 -45 45 -45 45 -45 45  0 0 -45 45 -45 45 -45 45  90];

for i = 1:length(theta)
    
    c = cosd(theta(i));
    s = sind(theta(i));
    
    Tinv{i} = [c^2                       s^2        -sqrt(2)*s*c;
                   s^2                       c^2         sqrt(2)*s*c;
                   sqrt(2)*s*c    -sqrt(2)*s*c    c^2 - s^2; ];
    
    T{i} = inv(Tinv{i});
    
    Q_bar{i} = Tinv{i}*(Q*T{i});

    Q11b(i) = Q_bar{i}(1,1);
    Q22b(i) = Q_bar{i}(2,2);
    Q12b(i) = Q_bar{i}(1,2);
    Q16b(i) = Q_bar{i}(1,3);
    Q26b(i) = Q_bar{i}(2,3);
    Q66b(i) = Q_bar{i}(3,3);
end

%% Vector h 
n = length(theta);
if mod(n,2) ~= 0
    n = n - 1;
end

% If all layers has the same thifness "e" and the ply is symetric
for i = 1:n/2
    h(i) = (n/2 - (i - 1))*e;
end

disp('The Vector of Quotas (h)')
h = [h 0 -flip(h)]

%% Matrix A, B and D

A11 = 0;  A12 = 0;  A22 = 0;  A66 = 0;
B11 = 0;  B12 = 0;  B22 = 0;  B66 = 0;
D11 = 0;  D12 = 0;  D22 = 0;  D66 = 0; D16 = 0; D26 = 0;

for k = 1:n
    
    % Matrix A
    A11 = Q11b(k)*abs(h(k)  - h(k+1)) + A11;
    A22 = Q22b(k)*abs(h(k)  - h(k+1)) + A22;
    A12 = Q12b(k)*abs(h(k)  - h(k+1)) + A12;
    A66 = Q66b(k)*abs(h(k)  - h(k+1)) + A66;

    % Matrix B
    B11 = Q11b(k)*(h(k)^2  - h(k+1)^2) + B11;
    B22 = Q22b(k)*(h(k)^2  - h(k+1)^2) + B22;
    B12 = Q12b(k)*(h(k)^2  - h(k+1)^2) + B12;
    B66 = Q66b(k)*(h(k)^2  - h(k+1)^2) + B66;

    % Matrix D
    D11 = Q11b(k)*(h(k)^3  - h(k+1)^3) + D11;
    D22 = Q22b(k)*(h(k)^3  - h(k+1)^3) + D22;
    D12 = Q12b(k)*(h(k)^3  - h(k+1)^3) + D12;
    D16 = Q16b(k)*(h(k)^3  - h(k+1)^3) + D16;
    D26 = Q26b(k)*(h(k)^3  - h(k+1)^3) + D26;
    D66 = Q66b(k)*(h(k)^3  - h(k+1)^3) + D66;
end

disp('The A matrix:')
A = [A11 A12 0;
       A12 A22 0;
       0 0 A66] 

disp('The B matrix:')
B = ([B11 B12 0; B12 B22 0; 0 0 B66])/2 

disp('The D matrix:')
D = ([D11 D12 D16; D12 D22 D26; D16 D26 D66])/3

%% ==================== Computing the minimum thickness =============================

disp('=========================')
disp('=====  FIRST METHOD ======')
disp('=========================')

h1 = e*n;

h_e = inv((1/h1)*A)*[Nx Ny sqrt(2)*Txy]';

for k = 1:n
    h_sigmaXY{k} = Q_bar{k}*h_e;
    
    h_sigmaLT{k} = T{k}*h_sigmaXY{k};
    
    %h_sigmaLT{k}(1)
    %h_sigmaLT{k}(2)
    %(h_sigmaLT{k}(3))/sqrt(2)
    
    hmin(k) = sqrt(((h_sigmaLT{k}(1))/(sigmaRL))^2 + ...
        ((h_sigmaLT{k}(2))/(sigmaRT))^2 - ...
        ((h_sigmaLT{k}(1))*(h_sigmaLT{k}(2)))/(sigmaRL^2) + ...
        (((h_sigmaLT{k}(3))/sqrt(2))/(tRLT))^2);
end

disp('First Method - The minimum thickness')
minimum_thickness = max(hmin)
message1 = strcat('First Method - The minimum thickness of laminate is', {' '}, num2str(minimum_thickness),{' '},'mm');
msgbox(message1)

disp('=========================')
disp('====  SECOND METHOD ====')
disp('=========================')

%%%%%%theta = [90 -45 45 -45 45 -45 45  0 0 -45 45 -45 45 -45 45  90];
hxA = Q_bar{1}*percentage_90 + Q_bar{2}*percentage_minus45 + Q_bar{3}*percentage_45 + Q_bar{8}*percentage_0;

h_e = inv(hxA)*[Nx Ny sqrt(2)*Txy]';

for k = 1:n
    h_sigmaXY{k} = Q_bar{k}*h_e;
    
    h_sigmaLT{k} = T{k}*h_sigmaXY{k};
    
    %h_sigmaLT{k}(1)
    %h_sigmaLT{k}(2)
    %(h_sigmaLT{k}(3))/sqrt(2)
    
    hmin(k) = sqrt(((h_sigmaLT{k}(1))/(sigmaRL))^2 + ...
        ((h_sigmaLT{k}(2))/(sigmaRT))^2 - ...
        ((h_sigmaLT{k}(1))*(h_sigmaLT{k}(2)))/(sigmaRL^2) + ...
        (((h_sigmaLT{k}(3))/sqrt(2))/(tRLT))^2);
end

disp('Second Method - The minimum thickness')
minimum_thickness = max(hmin)
message2 = strcat('Second Method - The minimum thickness of laminate is', {' '}, num2str(minimum_thickness),{' '},'mm');
msgbox(message2)

%% ==================== Plotting the Graph =============================
theta_graph = [0:1:90];

for i = 1:length(theta_graph)
    c = cosd(theta_graph(i));
    s = sind(theta_graph(i));
    
    Tinv{i} = [c^2                       s^2        -sqrt(2)*s*c;
                   s^2                       c^2         sqrt(2)*s*c;
                   sqrt(2)*s*c  -sqrt(2)*s*c    c^2 - s^2; ];
    
    T{i} = inv(Tinv{i});
    
    Q_rotate{i} = Tinv{i}*(Q*T{i});

    Qxx(i) = Q_rotate{i}(1,1);
    Qyy(i) = Q_rotate{i}(2,2);
    Qxy(i) = Q_rotate{i}(1,2);
end

angles = [90 45 0];
Qxx_plot = [Q11b(1) Q11b(3) Q11b(8)];
Qyy_plot = [Q22b(1) Q22b(3) Q22b(8)];
Qxy_plot = [Q12b(1) Q12b(3) Q12b(8)];


figure
fontsize = 16;
plot(theta_graph,Qxx,'b','LineWidth',2)
hold on
plot(theta_graph,Qyy,'k','LineWidth',2)
hold on
plot(theta_graph,Qxy,'r','LineWidth',2)
hold on
plot(angles,Qxx_plot,'m*','LineWidth',5)
hold on
plot(angles,Qyy_plot,'c*','LineWidth',5)
hold on
plot(angles,Qxy_plot,'g*','LineWidth',5)

legend('Qxx','Qyy','Qxy','Qxx for 0°,45° and 90°','Qyy for 0°,45° and 90°','Qxy for 0°,45° and 90°','FontSize',12)
xlabel('\theta [Degree]','FontSize',fontsize)
ylabel('Qs [MPa]','FontSize',fontsize)
title('Relation between Qxx, Qyy and Qxy and \theta','FontSize',fontsize)

grid on














