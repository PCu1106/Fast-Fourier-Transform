# Fast-Fourier-Transform

## complex 32-point FFT

(a)	theoretical derivations

<img src=./image/1.png width=60% />

Before implement 32-point FFT, I tried to program 8-point FFT. If I could find some regular rules for the counting, I can complete any N-point FFT. In the figure, I realized I could use three for loops to finish this work with three variables, stage, butterflies, and k. We could use these three variables to find the two indices we wanted to count and the value of coefficient W. Detail algorithm is in part(c).

(b)	flow diagram

<img src=./image/2.png width=60% />

We can change any N which $N=2^m$, m∈1,2,3,… For this question, we just assign N=32.

(c)	algorithms in Matlab

```matlab
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
```

(d)	complexity analyses

Each butterfly needs:
- 1 complex multiplication
- 2 complex additions

Each stage need: N/2 butterflies
- N/2   complex multiplications
- N complex additions

Total computation needs: $log{_2}{⁡N}$ stages
- $(Nlog{_2}{⁡N})/2$ complex multiplications
- $Nlog{_2}{⁡N}$ complex additions

## complex 64-point FFT

(a)	theoretical derivations

<img src=./image/4.png width=60% />

We can simply separate N-point DFT into two N/2-point DFT. Some are even numbers, and the others are odd numbers. (This is the basic idea of FFT.)

<img src=./image/5.png width=60% />
I drew the graph to help me implement the formula.

(b)	flow diagram

<img src=./image/6.png width=60% />

(c)	algorithms in Matlab

```matlab
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
```

(d)	complexity analyses

Each butterfly needs:
- 1 complex multiplication
- 2 complex additions

N/2-point FFT need: 
- $(N(log{_2⁡}{N}-1))/4$ complex multiplications
- $(N(log{_2}{⁡N}-1))/2$ complex additions

Total computation needs: 2 N/2-point FFTs & N/2 butterflies
- $(N(log{_2}{⁡N}-1))/2+N/2$ complex multiplications
- $N(log{_2}{⁡N}-1)+N$ complex additions

## double-real 64-point FFT

(a)	theoretical derivations

<img src=./image/8.png width=60% />

The main purpose is how to count two real number signals with only one FFT. We can use the symmetry properties in the figure. Let x, y be the two real number signals, z=x+jy, then count FFT of z. As a result, we can get:
$X[k]=FFT(x[n])=FFT(Re\lbrace z[n]\rbrace)=(Z[k]+Z^* [((-k))_N])/2$
$Y[k]=FFT(y[n])=FFT(Im\lbrace z[n]\rbrace)=(Z[k]-Z^* [((-k))_N])/2j$

(b)	flow diagram

<img src=./image/9.png width=60% />

(c)	algorithms in Matlab

```matlab
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
```

(d)	complexity analyses

creating z needs:
- N multiplications
- N additions

N-point FFT needs:
- $(N log{_2}{⁡N})/2$ complex multiplications
- $N log{_2}{⁡N}$ complex additions

Counting X needs:
- N complex multiplications
- N complex additions

Counting Y needs:
- N complex multiplications
- 2N complex additions

Total computation needs:
- $(N log{_2}{⁡N})/2+2N$ complex multiplications & N multiplications
- $N log{_2}{⁡N}+3N$ complex additions & N additions

## inverse complex 64-point FFT

<img src=./image/12.png width=60% />

We can simply use three step to implement inverse FFT.
Steps: FFT→/N→time reverse

(b)	flow diagram

<img src=./image/13.png width=60% />

(c)	algorithms in Matlab

```matlab
function x = IFFT64(X)
    N = 64;
    X1 = FFT64(X);
    x = X1./N;
    x(2:N) = fliplr(x(2:N));
end
```

(d)	complexity analyses

N-point FFT needs:
- $(N log{_2}{⁡N})/2$ complex multiplications
- $N log{_2}{⁡N}$ complex additions
/N needs:
- N complex multiplications
Total computation needs:
- $(N log{_2}{⁡N})/2+N$ complex multiplications
- $N log{_2}{⁡N}$ complex additions
