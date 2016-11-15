%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beginning of waveletTransformPlot.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function waveletTransformPlot(x, wavefileallName1, wavefileallName)
wname = 'db2';
Level = 4;

[C, L] = wavedec(x, Level, wname);
clf
subplot(1, 2, 1);
plot(x);
title('TV-curve')

cA3 = appcoef(C, L, wname, Level);

subplot(1, 2, 2);
plot(cA3);
title(['Approximation coefficients at level ', num2str(Level)])
print('-dsvg', wavefileallName1)

clf
for i = Level: - 1:1
    cD = detcoef(C, L, i);
    subplot(2, 2, Level - i + 1);
    plot(cD);
    title(['Detail coefficients at level ', num2str(i)])
end

print('-dsvg', wavefileallName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of waveletTransformPlot.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
