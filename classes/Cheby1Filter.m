classdef Cheby1Filter < IIRFilter
    %Cheby1Filter - Create Butterworth filter object.
    %   Fp , Fs can be scalars or 2dim vectors, in which case a band pass
    %   filter will be constructed.
    
    properties

    end
    
    methods
        function obj = Cheby1Filter(Fp , Fs , Ap , As , fs)
            if(nargin < 5)
                fs = [];
            end
            obj@IIRFilter(Fp , Fs , Ap , As , fs);
        end
    end
    
    methods (Access = protected)
        function [z,p,g] = designLowpassPrototype(obj)
            [OmegaProP , OmegaProS] = prototypeLowPassSpecs(obj);
            epsilon = sqrt(10^(obj.Ap/10)-1);
            N = Cheby1Filter.calculateOrder(OmegaProP, OmegaProS, obj.Ap , obj.As);
            Wc = OmegaProP;
            p = Cheby1Filter.calculatePoles(N , Wc ,epsilon);
            z = [];
            g = prod(-p);
            if(~mod(N , 2))
                g = g/sqrt(1+epsilon^2);
            end
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
        end
        
        function p = calculatePoles(N , Wc , epsilon)
            phi = 1/N*asinh(1/epsilon);
            theta = (2*(1:N)-1)*pi/(2*N) + pi/2;
            thetarev = wrev(theta);
            theta(1:2:N) = theta(1:ceil(N/2));
            theta(2:2:N) = thetarev(1:floor(N/2));
            p = Wc*(sinh(phi)*cos(theta) + 1i*cosh(phi)*sin(theta));
            p = transpose(p);
        end
        
    end
end

