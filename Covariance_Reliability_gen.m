clc;clear;close all
%%
% number of edges
N = 100;
NE = N*(N-1)/2;
% 
S = 0.2*ones(NE);
for i = 1:NE
    S(i,i) = 1;
end
save('data/Example_Covariance.mat','S');
writematrix(S,'data/Example_Covariance.csv');
%%
R_vec = 0.2 + (1-0.2)*rand(1,NE);
R = squareform(R_vec);

save('data/Example_Reliability.mat','R');
writematrix(R, 'data/Example_Reliability.csv');