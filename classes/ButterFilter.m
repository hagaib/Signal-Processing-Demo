classdef ButterFilter < IIRFilter
    %ButterFilter - Create Butterworth filter object.
    %   Fp , Fs can be scalars or 2dim vectors, in which case a band pass
    %   filter will be constructed.
    
    properties

    end
    
    methods
        function obj = ButterFilter(Fp , Fs , Ap , As , fs)
            if(nargin < 5)
                fs = [];
            end
            obj@IIRFilter(Fp , Fs , Ap , As , fs);
        end
    end
    
    methods (Access = protected)
        function [z,p,g] = designLowpassPrototype(obj)
            [OmegaProP , OmegaProS] = prototypeLowPassSpecs(obj);
            N = ButterFilter.calculateOrder(OmegaProP, OmegaProS, obj.Ap , obj.As);
            Wc = ButterFilter.calculate3dbCutoff(N , OmegaProS , obj.As);
            p = ButterFilter.calculatePoles(N , Wc);
            z = [];
            g = Wc^N;
%             b = g*poly(z);
%             a = poly(p);
%             [h,w] = freqs(b , a , 4096);
%             figure(1)
%             plot(w , 20*log10(abs(h)));
        end
    end
    
    
    methods (Static)
        
        function Wc = calculate3dbCutoff(N , Ws , As)
            Wc = Ws*(10^(As/10)-1)^(-1/(2*N));
        end
        
        function N = calculateOrder(Wp , Ws , Ap , As)
            N = ceil(0.5 * (log10(10^(As/10)-1) - log10(10^(Ap/10)-1) ) / log10(Ws/Wp));
        end
        
        function p = calculatePoles(N , Wc)
            theta = (2*(1:N)-1)*pi/(2*N) + pi/2;
            thetarev = wrev(theta);
            theta(1:2:N) = theta(1:ceil(N/2));
            theta(2:2:N) = thetarev(1:floor(N/2));
            p = Wc*exp(1i*theta);
            p = transpose(p);
        end
        
    end
end

