clear

% T = readtable('143fires/abneydf.csv');
% 
% a = table2array(T);
% 
% t = a(:,2)/a(end,2);
% A = a(:,1)/a(end,1);
% 
% plot(t,A,'LineWidth',2)


P = '143fires'; 
S = dir(fullfile(P,'*.csv')); 
N = numel(S); %143 fires
for k = 1:N
    F = fullfile(P,S(k).name);
    t = readtable(F);
    h = height(t);
    H(k) = h;
    M = sum(H>0);
end
% max(H) = 89;114

T = nan(89,2,M);
Names = strings(M,1);
i = 1;
for k = 1:N
    F = fullfile(P,S(k).name);    
    t = readtable(F);
    a = table2array(t);
    n = size(a,1);
    if n > 0
        name = S(k).name;
        Names(k) = name(1:end-4);
        area_dim = a(:,1)/a(end,1);
        time_dim = a(:,2)/a(end,2);
        T(1:n,1,i) = time_dim;
        T(1:n,2,i) = area_dim;
        i = i+1;
    end
    clearvars('time_dim'); clearvars('area_dim')
end

save("selectedfires.mat",'T','Names')

