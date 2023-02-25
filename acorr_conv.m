% #2-1
x = zeros(1, 64);
x(1:32) = [3:3:96];
h = zeros(1, 64);
h(1:32) = [30:-2:-32];
conv1 = direct_conv(x, h);
figure(1);
subplot(2,1,1);
stem(linspace(1,63,63), conv1);
title('a direct real convolution program');
conv2 = DRFFT64_conv(x, h);
subplot(2,1,2);
stem(linspace(1,63,63), conv2);
title('a real convolution program by calling DRFFT64(x, y) once');
error = conv1 - conv2;
disp(error);

% #2-2
x = zeros(1, 64);
n = [1:32];
x(1:32) = n.*(-1).^n;
acorr1 = direct_acorr(x);
figure(2);
subplot(2,1,1);
stem(linspace(1,63,63), acorr1);
title('a direct real autocorrection computation program');
acorr2 = DRFFT64_acorr(x);
subplot(2,1,2);
stem(linspace(1,63,63), acorr2);
title('a computation program by calling DRFFT64 (x, y) once');
error = acorr1 - acorr2;
disp(error);

function acorr = DRFFT64_acorr(x)
    x2 = zeros(1, 64);
    x2(33:64) = x(1:32);
    [X1, X2] = DRFFT64(x, x2);
    ACORR = X1.*conj(X2);
    ans = IFFT64(ACORR);
    acorr = ans(2:64);
end

function acorr = direct_acorr(x)
    x2 = zeros(1, 64);
    x2(33:64) = x(1:32);
    ans = zeros(1, 2*length(x)-1);
    for n = 1:64
        for m = 1:n
            ans(n) = ans(n) + (x(m) * x2(end-n+m));
        end
    end
    acorr = ans(1:63);
end

function conv = DRFFT64_conv(x, h)
    [X, H] = DRFFT64(x, h);
    CONV = X.*H;
    ans = IFFT64(CONV);
    conv = ans(1:63);
end

function conv = direct_conv(x, h)
    h_rev = fliplr(h);
    ans = zeros(1, length(x)+length(h)-1);
    for n = 1:64
        for m = 1:n
            ans(n) = ans(n) + (x(m) * h_rev(end-n+m));
        end
    end
    conv = ans(1:63);
end

function x = IFFT64(X)
    N = 64;
    X1 = FFT64(X);
    x = X1./N;
    x(2:N) = fliplr(x(2:N));
end

function [X, Y] = DRFFT64(x, y)
    N = 64;
    z = x+1i*y;
    Z = FFT64(z);
    Zflip(1) = Z(1);
    Zflip(2:N) = fliplr(Z(2:N));
    ZZ = conj(Zflip);
    X = (Z+ZZ)/2;
    Y = -1i*(Z-ZZ)/2;
end

function X = FFT64(x)
    N = 64;
    x1 = x(1:2:N);
    x2 = x(2:2:N);
    G = FFT32(x1);
    H = FFT32(x2);
    X = zeros(1,N);
    for k = 0:(N/2)-1
       temp1 = G(k+1);
       temp2 = H(k+1)*W(k,N);
       X(k+1) = temp1+temp2;
       X(k+1+(N/2)) = temp1-temp2;
    end
end

function X = FFT32(x)
    N = 32;
    indices = 0 : length(x)-1;
    new_indices = bin2dec(fliplr(dec2bin(indices, log2(N))));
    X = x(new_indices + 1);
    for stage = 0:log2(N)-1
        for butterflies = 0:(N/2^(stage+1))-1
            for k = 0:(2^stage)-1
                first = butterflies*(2^(stage+1))+k+1;
                second = first+(2^stage);
                temp1 = X(first);
                temp2 = X(second)*W(k,2^(stage+1));
                X(first) = temp1+temp2;
                X(second) = temp1-temp2;
            end
        end
    end
end

function Wvalue = W(k,N)
    Wvalue = exp(-1i * k * 2 * pi / N);
end