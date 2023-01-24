
%%% Arian Velayati, PhD
%%%% This script is used to find the Young's modulus of rocks using the 50%
%%%% max peak stress method

clc; clear; close;
%% Input

L_o = 1.9990; % Sample length (in)
D_o = 0.99825; % sample diameter (in)
ea2_M = load('ea2_M.txt'); % Displacement (mil)
Sd_M = load('Sd_M.txt'); % Deviatoric Stress (psi)
et_M = load('et_M.txt'); % Radial strain (mil)

%% Calculations

ea2_M = ea2_M./L_o; % Displacement to ENG strain (mil/in)
ea2_M = ea2_M*0.001; % (mil/mil)
Sd_m = max(Sd_M); % Max Sd value (psi)
Sd_mh = 0.5*Sd_m; % Half Sd_m (psi)
et_M = et_M./D_o; % radial strain (mil/in)
et_M = et_M.*0.001; % radial strain (mil/mil)

r = find(Sd_M>Sd_mh); % Find the row in which the 0.5 half peak stress can be found

% YM
n = 3; % arbitrary row distancing from the 50%
f = polyfit(ea2_M(r(1)-n:(r(1)+n)), Sd_M(r(1)-n:r(1)+n),1); %Fitting a linear tanget around 50% of peak s1
fout = f;
for i = 1:10
    % f = polyfit(ea((r(1)+n):r(1)), Sd((r(1)+n):r(1)),1)
    n = n + 180;
f = polyfit(ea2_M(r(1)-n:(r(1)+n)), Sd_M(r(1)-n:r(1)+n),1);
fout = [fout;f];
if (fout(i+1,1)-fout(i,1))/fout(i+1,1) < 0.01
    disp('Optimal n is: '); disp(n); 
    disp('E (GPa) is: '); disp(f(:,1)/(145*1000))
    E = f(:,1)/(145*1000); %GPa
     break
end
end

figure(1)
scatter(ea2_M,Sd_M)
hold on 
x = linspace(ea2_M(r(1)-n),ea2_M(r(1)+n),100);
y = f(1).*x+f(2);
plot(x,y,'r','LineWidth',5)
xlabel('Eng Strain (mil/mil)')
ylabel('Stress (psi)')
legend('Experimental Data','50% tangent method')

% PR
n = 10; % arbitrary row distancing from the 
f = polyfit(ea2_M(r(1)-n:(r(1)+n)), et_M(r(1)-n:r(1)+n),1); %Fitting a linear tanget around 50% of peak s1
fout = f;
for i = 1:10
    % f = polyfit(ea((r(1)+n):r(1)), Sd((r(1)+n):r(1)),1)
    n = n + 50;
f = polyfit(ea2_M(r(1)-n:(r(1)+n)), et_M(r(1)-n:r(1)+n),1);
fout = [fout;f];
if (fout(i+1,1)-fout(i,1))/fout(i+1,1) < 0.01
    disp('Optimal n is: '); disp(n); 
    disp('Poisson ratio: '); disp(-1*f(end,1))
     break
end
end

%% plot
 
figure(2)
scatter(ea2_M,et_M)
hold on 
x = linspace(ea2_M(r(1)-n),ea2_M(r(1)+n),100);
y = f(1).*x+f(2);
plot(x,y,'r','LineWidth',5)
xlabel('Eng Axial Strain (mil/mil)')
ylabel('Eng Lateral Strain (mil/mil)')
legend('Experimental Data','50% tangent method')

figure(3)
scatter(ea2_M,Sd_M)
hold on
scatter(et_M,Sd_M)

%% Method2 to find PR
n = 1; % arbitrary row distancing from the 
f = polyfit(et_M(r(1)-n:r(1)+n), Sd_M(r(1)-n:r(1)+n),1); %Fitting a linear tanget around 50% of peak s1
fout = f;
for i = 1:10
    % f = polyfit(ea((r(1)+n):r(1)), Sd((r(1)+n):r(1)),1)
    n = n + 50;
f = polyfit(et_M(r(1)-n:r(1)+n), Sd_M(r(1)-n:r(1)+n),1);
fout = [fout;f];
if (fout(i+1,1)-fout(i,1))/fout(i+1,1) < 0.01
    disp('Optimal n is: '); disp(n); 
    disp('Poisson ratio2: '); disp((E*1000*145)/(-1*f(end,1)))
         break
end
end
v = (E*1000*145)/(-1*f(end,1)); % PR2

figure(4)
scatter(et_M,Sd_M)
hold on
x = linspace(et_M(r(1)-n),et_M(r(1)+n),100);
y = f(1).*x+f(2);
plot(x,y,'r','LineWidth',5)
xlabel('Eng lateral Strain (mil/mil)')
ylabel('Deviatoric Stress (psi)')
legend('Experimental Data','50% tangent method')

%% Elastic Constants

k = E/(3*(1-2*v)); % GPa
G = E/(2*(1+v)); % GPa

  T = table(E,G,v,k,'VariableNames',{'YM_GPa','SM_GPa','PR','BM_GPa'})
  
    writetable(T,'ULT_Moduli.csv')