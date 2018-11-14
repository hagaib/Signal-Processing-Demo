classdef Cheby2Filter < IIRFilter
    %Cheby2Filter - Create Chebyshev 2 filter object.
    %   Fp , Fs can be scalars or 2dim vectors, in which case a band pass
    %   filter will be constructed.
    
    properties

    end
    
    methods
        function obj = Cheby2Filter(Fp , Fs , Ap , As , fs)
            if(nargin < 5)
                fs = [];
            end
            obj@IIRFilter(Fp , Fs , Ap , As , fs);
        end
    end
    
    methods (Access = protected)
        function [z,p,g] = designLowpassPrototype(obj)
            [OmegaProP , OmegaProS] = prototypeLowPassSpecs(obj);
            A = 10^(obj.As/20);
            N = Cheby2Filter.calculateOrder(OmegaProP, OmegaProS, obj.Ap , obj.As);
            Wc = OmegaProS;
            p = Cheby2Filter.calculatePoles(N , Wc , A);
            z = Cheby2Filter.calculateZeros(N , Wc);
            g = real(prod(-p)/prod(-z));
%             b = g*poly(z);
%             a = poly(p);
%             [h,w] = freqs(b , a , 4096);
%             figure(1)
%             plot(w , 20*log10(abs(h)));
        end
    end
    
    
    methods (Static)
        
        function N = calculateOrder(Wp , Ws , Ap , As)
            alpha = Ws/Wp;
            beta = sqrt( ( 10^(As/10)-1 ) / ( 10^(Ap/10)-1) );
            N = ceil(acosh(beta)/acosh(alpha));
            if(~mod(N,2)) , N=N+1; end
        end
        
        function z = calculateZeros(N , Wc)
            theta = (2*(1:N)-1)/(2*N);
            thetarev = wrev(theta);
            theta(1:2:N) = theta(1:ceil(N/2));
            theta(2:2:N) = thetarev(1:floor(N/2));
            % Remove zero at infinity
            theta = theta(theta ~= 0.5);
            theta = pi * theta;
            z = 1i*Wc./cos(theta);
            z = transpose(z);
        end
        
        function p = calculatePoles(N , Wc , A)
            cheby1poles = Cheby1Filter.calculatePoles(N , Wc , 1/sqrt(A^2-1) );
            p = Wc^2./cheby1poles;
        end
        
    end
end

