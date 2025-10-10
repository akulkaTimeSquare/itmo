temp = imread("log_edited - Copy.png");
img = double(temp)/255;
fourier = fftshift(fft2(img));
ab = abs(fourier); % абсолютные значения
an = angle(fourier); % фазы
l = log(ab+1);
m = max(l(:));
ul = l/m; % нормировка логарифма
imwrite(ul, "edited_1.png");
imshow(ul, []);