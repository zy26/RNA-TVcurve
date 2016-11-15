% x=rand(200,1);
% y=rand(238,1);

function s = waveletSimilarity(x, y, wname, Level)
% %% x,y 原始信号
% %%wname: 选择的小波  %%默认是db2 静态小波变换
% %%level：小波变换的层数: 默认是Maximum wavelet decomposition level
%
%
%

if nargin == 2,
    wname = 'db2';
    Timex = length(x);
    Timey = length(y);
    Levelx = wmaxlev(Timex, wname);
    Levely = wmaxlev(Timey, wname);
    Level = min(Levelx, Levely);
    disp(['Maximum wavelet decomposition level:  ', num2str(Level)]);
end
if nargin == 3,
    Timex = length(x);
    Timey = length(y);
    Levelx = wmaxlev(Timex, wname);
    Levely = wmaxlev(Timey, wname);
    Level = min(Levelx, Levely);
    disp(['Maximum wavelet decomposition level:  ', num2str(Level)]);
end

if nargin == 4,
    disp(['Wavelet decomposition level:  ', num2str(Level)]);
    Timex = length(x);
    Timey = length(y);
    wmaxlev(Timex, wname);
    wmaxlev(Timey, wname);
end
disp(['Wavelet name:  ', wname]);

%   Solve suggested length
pow = 2 ^ Level;
if rem(Timex, pow) > 0
    sugLengthx = ceil(Timex / pow) * pow;
else
    sugLengthx = Timex;
end
if rem(Timey, pow) > 0
    sugLengthy = ceil(Timey / pow) * pow;
else
    sugLengthy = Timey;
end

disp(['Note: 2^Level has to divide the length of the signal']);
disp(['Suggested length of  x: ' int2str(sugLengthx), '. Suggested length of y :', int2str(sugLengthy)]);

if sugLengthy > sugLengthx
    Lenx = sugLengthy - Timex;
    Leny = sugLengthy - Timey;
else
    Lenx = sugLengthx - Timex;
    Leny = sugLengthx - Timey;
end
%    Extend vector or matrix

x1 = wextend('1', 'zpd', x, Lenx, 'r');
y1 = wextend('1', 'zpd', y, Leny, 'r');

[swtlowcoeX, swthighcoeX] = swt(x1, Level, wname);
[swtlowcoeY, swthighcoeY] = swt(y1, Level, wname);

lowcorr = abs(corr(swtlowcoeX(Level, :)',swtlowcoeY(Level,:)'));

weight(Level + 1) = 2 ^ ((Level - 1) / 2);
s = lowcorr * weight(Level + 1);
highcorr = zeros(Level);
for m = 1:Level
    highcorr(m) = abs(corr(swthighcoeX(m, :)',swthighcoeY(m,:)'));
    weight(m) = 2 ^ ((m - 1) / 2);
    s = s + highcorr(m) * weight(m);
end
s = s / sum(weight);