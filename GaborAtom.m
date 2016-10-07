function [g,parameters,num] = GaborAtom(N)
num = 0;
 M = 12*(4*N*log2(N)+12*N+log2(N)-3);
 g = zeros(M,N);
 parameters = zeros(M,4);
for j = 0:log2(N)
     s = 2^j;
    for p = 0:N*2^(-j+1)
         u = s*p*1/2;
        for k = 0:2^(j+1)
            v = 2^(-j)*pi*k;
            for i = 0:11
                w = i*pi/6;
                num = num +1;
                g(num,:) = getGaborAtom(N,s,u,v,w);
                parameters(num,:) = [s u v w];
            end  
        end
    end
end


%%-----------------------------------
function atom = getGaborAtom(N,scale,timeShift,frequency,phase)

atom =zeros(N,1);
for n=1:N
    atom(n,1) = (1/sqrt(scale))*exp(-pi*(n-timeShift)^2/scale^2) * cos(frequency*n+phase);
end
atom = (1/norm(atom)) .* atom;   %Normalization
