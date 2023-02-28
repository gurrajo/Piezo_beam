clear all; close all;clc
a = linspace(0,10000,100);
b = linspace(0,10000,100);
[A,B] = meshgrid(a,b);
odds_1 = 1.6;
odds_2 = 2.70;

f_1 = (odds_1-1)*A - B;
f_2 = (odds_2-1)*B - A;
C = min(f_1,f_2);
figure
surf(A,B,C)

[M,ind1] = max(C);
[max_val,ind2] = max(M);
bet_sum = A(ind1(ind2),ind2) +  B(ind1(ind2),ind2);
gain = max_val/bet_sum;

disp("Kronor på lag A " + A(ind1(ind2),ind2))
disp("Kronor på lag B " + B(ind1(ind2),ind2))
disp("Vinst om lag A vinner: " + f_1(ind1(ind2), ind2));
disp("Vinst om lag B vinner: " + f_2(ind1(ind2), ind2));
disp("Procentuel vinst: " + gain*100 + " %")
