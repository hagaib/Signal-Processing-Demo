classdef FIRFilter < Filter
    %FIRFILTER generates FIR filters using the windowing method to filter Signal objects.
    % only Hamming and Hann windows are supported 
    % only even order filters are supported
    
    properties
        window  % window name (string)
        winid   % id of window according to WindowGenerator class
        Wc      % cutoff frequency (normalized)
    end
    
    methods
        function obj = FIRFilter(window , Fp , Fs , Ap , As , fs)  
            
            obj@Filter(Fp , Fs);
            
            if(nargin <6)
                fs = 2;
            end
            obj.Fp = Fp;
            obj.Fs = Fs;
            obj.Wp = Fp/fs*2;
            obj.Ws = Fs/fs*2;
            obj.Ap = Ap;
            obj.As = As;
            
            if(isequal(window , 'hamming'))
                obj.window = 'hamming';
                obj.winid = 7;
            elseif(isequal(window , 'hanning'))
                obj.window = 'hanning';
                obj.winid = 4;
            end
            
            delta_w = max(abs(obj.Ws-obj.Wp));
            obj.Wc = (obj.Wp+obj.Ws)/2;
            obj.N = filterorder(obj , delta_w);
            createFilter(obj);
        end
        
        
        function plot(obj , fignum)
            [h,w] = freqz(obj.b , obj.a , 4096);
            figure(fignum)
            subplot(3,1,1)
            plot(w/(pi) , 10*log10(h.*conj(h)))
            title([class(obj) , ' magnitude response'])
            subplot(3,1,2)
            plot(w/(pi) , unwrap(angle(h)))
            title('phase response')
            subplot(3,1,3)
            plot(roots(obj.b) , 'o');
            hold on
            plot(cos(0:0.1:2*pi) , sin(0:0.1:2*pi));
            hold off
            title('zeros and poles on the Z plane');
            legend('zeros')
        end
        
        function n = filterorder(obj , delta_w)
            windowtype = obj.window;
            if(isequal(windowtype,'hamming'))
                n = ceil(8/delta_w);
            elseif (isequal(windowtype,'hanning'))
                n = ceil(8/delta_w);
            elseif(isequal(windowtype,'rectwin'))
                n = ceil(4/delta_w);
            end
            
            if(length(obj.Ws)>1)
                n=floor(n*1.05);
            end
            if(mod(n,2)==1) , n=n+1; end
        end
        
        function createFilter(obj)
            iideal = ideal_impulse_response(obj);
            wingen = WindowGenerator(obj.winid , obj.N+1);
            w = wingen.window';
            obj.b = iideal .* w;
            obj.a = 1;
        end
        
        function iideal = ideal_impulse_response(obj)
            W = obj.Wc;
            N = obj.N;
            if(length(obj.Wc)>1) , W = obj.Wc(2);end;
            
            t = -N/2:N/2;
            iideal = sin(pi*W.*t)./(pi*t);
            iideal(N/2+1)=W;
            
            if(length(obj.Wc)>1)
                %bandpass
                W = obj.Wc(1);
                iideal2 = sin(pi*W.*t)./(pi*t);
                iideal2(N/2+1)=W;
                
                iideal = iideal - iideal2;
                
            elseif(obj.Wp(1)>obj.Ws(1))
                %highpass
                iideal = -iideal;
                iideal(N/2+1)=iideal(N/2+1)+1;
            end
        end
        
    end
    
    
end

