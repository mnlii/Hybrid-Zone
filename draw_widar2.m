close all;
letters = {I, N, F4, O6, C2, O9, M15};
for i = 1:length(letters)
    tmp = letters{i};
    figure;
    plot(tmp(1,:), tmp(2,:));
end