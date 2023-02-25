x = [1 2 3 4 5 6 7 8 9 10 11 0 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 14 14 15 16 1 2 3 4 5 6 7 8 7 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
y = [9 2 3 4 5 6 7 8 9 18 11 0 13 14 15 16 1 2 3 4 5 6 7 81 9 10 11 20 14 14 15 16 1 2 3 4 5 6 7 8 7 10 11 12 15 14 15 16 1 2 3 4 5 6 27 8 9 10 11 12 13 14 15 16];
[X1, Y1] = DRFFT64(x, y);
X2 = FFT64(x);
Y2 = FFT64(y);
error1 = X1 - X2;
error2 = Y1 - Y2;
disp(error1);
disp(error2);
%{
Y = [1 1 2 2 3i 3+4i 4 4 5i 5 6 6 7+7i 7 8i 8 1i 1 2-9i 2 3i 3 4 4 5i 5 6 6 7 7 8i 8 1 1 2 2 3i 3+4i 4 4 5i 5 6i 6 7+7i 7 8i 8 1i 1 2-8i 2 3i 3 4 4 5i 5 6 6 7+7i 7 8i 8];
y = ifft(Y,64);
disp(y);
y = IFFT64(Y);
disp(y);
%}

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