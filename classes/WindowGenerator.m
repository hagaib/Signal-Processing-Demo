classdef WindowGenerator < handle
    % only windows with uneven number of samples are supported
    % Generates windows according to name or number:
    %1 Rect
    %2 Triangle
    %3 Hanning 1
    %4 Hanning 2
    %5 Hanning 3
    %7 Hamming
    %8 Riesz
    %9 Riemann
    %10 de la Valle Poussin
    %11 Tukey 0.25
    %12 Tukey 0.50
    %13 Tukey 0.75
    %14 Bohman
    %15 Poisson a=2
    %16 Poisson a=3
    %18 Cauchy 3
    %19 Cauchy 4
    %21 Gaussian a=2.5
    %22 Gaussian a=3.0
    %23 Gaussian a=3.5
    %24 Kaiser-Bessel a=2.0
    %25 Kaiser-Bessel a=2.5
    %26 Kaiser-Bessel a=3.0
    %28 Exact Blackman
    %29 Blackman
    properties
        window      % samples of window function
        name
        id
        N           % window sample count
        x           % window support
    end
    
    properties (Constant)
        winCount = 29
    end
    
    methods (Access = public)
        function obj = WindowGenerator(varargin)
            if(nargin<1)
                obj.name = nan;
                obj.window = nan;
                return;
            elseif(nargin<3)
                type = 'noteven';
            else
                type = varargin{3};
            end
            generate(obj , varargin{1} ,varargin{2} , type);
        end
        
        function list = allWindowNames(obj)
            list = cell(obj.winCount , 1);
            for i=1:obj.winCount
                generate(obj , i , 3);
                list{i} = obj.name;
            end
        end
        
        function win = generate(obj , winnum , N , type)
            obj.id = winnum;
            iseven = ~mod(N,2);
            if(iseven) , N=N+1; end
            obj.N = N;
            x = floor(-(N-1)/2:(N-1)/2);
            x=x';
            
            switch winnum
                case 1
                    obj.name = 'Rectangle';
                    winfunc = @(n)ones(length(n),1);
                case 2
                    obj.name = 'Triangle';
                    winfunc = @(n) 1 - abs(n)/((N-1)/2);
                case 3
                    obj.name = 'Hanning 1 (Cosine^1)';
                    winfunc = @(n) cos(n/N*pi);
                case 4
                    obj.name = 'Hanning 2 (Real Hanning)';
                    winfunc = @(n) cos(n/N*pi).^2;
                case 5
                    obj.name = 'Hanning 3 (Cosine^3)';
                    winfunc = @(n) cos(n/N*pi).^3;
                case 7
                    obj.name = 'Hamming';
                    alpha = 0.54;
                    winfunc = @(n) alpha + (1-alpha)*cos(2*pi/N*n);
                case 8
                    obj.name = 'Riesz';
                    winfunc = @(n) 1 - abs(n/(N-1)*2).^2;
                case 9
                    obj.name = 'Riemann';
                    winfunc = @(n) sinc(n/(N-1)*2);
                case 10
                    obj.name = 'de  la  Valle  Poussin';
                    winfunc = @(n) delavallewin(obj , n);
                case 11
                    obj.name = 'Tukey \alpha = 0.25';
                    winfunc = @(n) tukeywin(obj , n , 0.25);
                case 12
                    obj.name = 'Tukey \alpha = 0.5';
                    winfunc = @(n) tukeywin(obj , n , 0.5);
                case 13
                    obj.name = 'Tukey \alpha = 0.75';
                    winfunc = @(n) tukeywin(obj , n , 0.75);
                case 14
                    obj.name = 'Bohman';
                    winfunc = @(n) ( 1 - abs(n)./(N-1)*2).*cos(pi.*abs(n)./(N-1)*2) + 1/pi*sin(pi.*abs(n)./(N-1)*2);
                case 15
                    obj.name = 'Poissond \alpha = 2';
                    winfunc = @(n) poissonwin(obj , n , 2);
                case 16
                    obj.name = 'Poisson \alpha = 3';
                    winfunc = @(n) poissonwin(obj , n , 3);
                case 18
                    obj.name = 'Cauchy \alpha = 3';
                    winfunc = @(n) cauchywin(obj , n , 3);
                case 19
                    obj.name = 'Cauchy \alpha = 4';
                    winfunc = @(n) cauchywin(obj , n , 4);
                case 21
                    obj.name = 'Gauss \alpha = 2.5';
                    winfunc = @(n) gaussianwin(obj , n , 2.5);
                case 22
                    obj.name = 'Gauss \alpha = 3.0';
                    winfunc = @(n) gaussianwin(obj , n , 3.0);
                case 23
                    obj.name = 'Gauss \alpha = 3.5';
                    winfunc = @(n) gaussianwin(obj , n , 3.5);
                case 24
                    obj.name = 'Kaiser-Bessel \alpha = 2.0';
                    winfunc = @(n) kaiserbessel(obj , n , 2.0);
                case 25
                    obj.name = 'Kaiser-Bessel \alpha = 2.5';
                    winfunc = @(n) kaiserbessel(obj , n , 2.5);
                case 26
                    obj.name = 'Kaiser-Bessel \alpha = 3.0';
                    winfunc = @(n) kaiserbessel(obj , n , 3.0);                    
                case 28
                    obj.name = 'Exact Blackman';
                    a0 = 0.42659071;
                    a1 = 0.49656062;
                    a2 = 0.07684867;
                    winfunc = @(n) a0 + a1*cos(2*pi/N.*n) + a2*cos(2*pi/N*2.*n);
                case 29
                    obj.name = 'Blackman';
                    a0 = 0.42;
                    a1 = 0.5;
                    a2 = 0.08;
                    winfunc = @(n) a0 + a1*cos(2*pi/N.*n) + a2*cos(2*pi/N*2.*n);
                otherwise
                    disp(['Window number ' , winnum , ' is undefined.']);
                    obj.name = nan;
                    obj.window = nan;
                    win = nan;
                    return
                    
                    
            end
            
            win = winfunc(x);
            
            if(nargin>3 && strcmp(type, 'dfteven'))
                % DFT - Even Window:
                win = win(1:end-1);
            end
            
            obj.x = x;
            obj.window = win;
        end
        
        function plot(obj, fignum)
            if(~isnan(obj.window))
                if(nargin<2)
                    figure
                else
                    figure(fignum)
                end
                N = obj.N;
                
                subplot(2,1,1)
                plot(obj.x(1:length(obj.window)) , obj.window)
                
                title([obj.name , ' Window']);
                
                subplot(2,1,2)
                fftN = 2^nextpow2(obj.N*100);
                spectwin = abs(fft(obj.window, fftN));
                A = spectwin(1);
                plot((0:(fftN-1))/fftN*2*pi-pi , 20*log10(fftshift(spectwin)./A))
                
                title('DTFT (interpolated)');
            end
        end
    end
    
    methods(Access = private)
        function y = delavallewin(obj , n)
            N = obj.N;
            
            y = (abs(n)<=N/4) .* (1 - 6.*(n./N.*2).^2 .* (1 - abs(n)./N.*2)) ...
                + ...
                (abs(n)>N/4) .* (2.*(1 - abs(n)./N.*2).^3) ; 
        end
        
        function y = poissonwin(obj , n ,alpha)
            N = obj.N;
            
            y = exp(-alpha .* abs(n) / (N-1)*2);
        end
        
        function y = kaiserbessel(obj , n , alpha)
            N = obj.N;
            b1 = besseli(0 , pi*alpha *sqrt(1 - (n/(N-1)*2).^2) );
            b2 = besseli(0 , pi*alpha);
            y = b1 ./ b2;
        end
        
        function y = cauchywin(obj , n ,alpha)
            N = obj.N;
            
            y = 1 ./ (1 + (alpha .* n ./(N-1)*2).^2 );
        end
        
        function y = gaussianwin(obj , n ,alpha)
            N = obj.N;
            
            y = exp(-0.5 * (alpha .* n ./ (N-1)*2).^2);
        end
        
        function y = tukeywin(obj , n , alpha)
            N = obj.N;
            alpha = 1-alpha;
            
            
            y1 = (abs(n)<=N/2*alpha);
            y2 = (abs(n)>N/2*alpha) .* 0.5 .* (1 + cos(pi .* (abs(n) - alpha.*(N/2)) ./ ((1-alpha)*(N/2)) ) ); 
            
            if(alpha == 1)
                y2 = zeros(size(y2));
            end
            
            y = y1+y2;
        end
        
    end
    
end

