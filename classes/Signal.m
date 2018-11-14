classdef Signal < matlab.mixin.Copyable
    
    % Init with c'tor, supplying an 8 length array
    % every signal is made of 
    
    % %Example:
    %s = Signal([2 ,2.3 , 0.4  , 0 , 0 , 0 , 0]);
    % % Every signal is made of several components. They are combined
    % additively
    %s.addcomp([3 ,17.7 , 1.2 , 2.5 , 9.5 , 8 , 2]);
    % %Generate by supplying sample rate and time:
    % fs = 44.1; %khz
    % time = 3000; %ms
    %s.generate(fs,time);
    
    properties
        title
        components % cell array
        xx
        fs
    end
    
    methods
        
        %% Constructor
        function obj = Signal(arr)
            if(nargin<1)
                obj.components = [];
            else
                if(ischar(arr))
                    obj.readCompsFromFile(arr);
                else
                    obj.components{1} = Component(arr);
                end
            end
        end
        
        %% Initialization
        
        function addComp(obj, arr)
            l = length(obj.components);
            obj.components{l+1} = Component(arr);
        end
        
        function generate(obj,sample_rate , time)
            obj.xx = zeros(size(floor(time*sample_rate),1));
            for i=1:length(obj.components)
                x = obj.components{i}.generate(sample_rate , time);
                obj.xx = obj.xx + x;
            end
            obj.fs = sample_rate;
        end
        
        
        function timevect = timevector(obj)
            if(isempty(obj.xx))
                timevect = [];
            else
                timevect = (0:(length(obj.xx)-1))/obj.fs';
            end
            
        end
        
        
        %% Plot
        function plot(obj , fignum , clist)
            if(nargin < 3)
                clist = 1:length(obj.components);
            end
            if(nargin<2)
                figure
            else
                figure(fignum);
            end
            t= obj.timevector;
            subplot(2,1,1);
            plot(1,1 ,t , obj.xx);
            hold on
            title(obj.title)
            for i=1:length(clist)
                if(length(obj.components{clist(i)}.xx) == length(t))
                    plot(t , obj.components{clist(i)}.xx);
                end
            end
            hold off
%             legend(lstrings , obj.title);
            subplot(2,1,2);
            plot(1,2, (0:(length(t)-1))/(length(t)-1)*obj.fs, abs(fft(obj.xx)))
            
%             legend(lstrings , obj.title);
        end
        
        %% I/O Methods
        function readCompsFromFile(obj , filename)
            input_data = csvread(filename);
            comps = size(input_data,2);
            for j=1:comps
                obj.addComp(input_data(:,j));
            end
        end
        
        function writeToFile(obj,filename)
            signal = obj.xx;
            fs = obj.fs;
            save(filename , 'signal' , 'fs' , 'obj');
            disp(['Writing to File: ' ,obj.title]);
        end
        
        %% Sampling and Complex Phase
        
        function resample(obj , factor)
            obj.xx = obj.xx(1:factor:end);
            obj.fs = obj.fs/factor;
            obj.title = ['sampled ' , obj.title];
        end
        
        function iq_env = IQ(obj , midfreq , bandwidth)
            %fpd = center frequency for demodulation
            t = obj.timevector;
            lpf = ButterFilter(bandwidth , bandwidth+0.4 , 0.4 , 34.7 , obj.fs);

            Imix = cos(2*pi*midfreq*t)' .* obj.xx;
            Qmix = sin(2*pi*midfreq*t)' .* obj.xx;

            I = lpf.filter3(Imix);
            Q = lpf.filter3(Qmix);

            IQdem = I - 1i*Q;
            iq_env = copy(obj);
            iq_env.xx = IQdem;
            iq_env.components = [];
            iq_env.title = ['IQ demodulated ' , iq_env.title];
        end
        
        %% Spectral Analysis
        function [t , f , s] = spectralAnalysis(obj , winnum, windur , zeropad, overlap, fignum)
            winlen = floor(windur*obj.fs);
            wingen = WindowGenerator();
            winlen = winlen + mod(winlen,2);
            window = wingen.generate(winnum , winlen+1 , 'dfteven');
            dftn = length(window);
            if(nargin >= 4 && ~isempty(zeropad))
                dftn = zeropad;
            end
            if(nargin < 5 || isempty(overlap))
                overlap = 0;
            end
            step = floor((1-overlap)*winlen);
            disp(['Spectral Analysis using a ', num2str(windur) , ' ms ' , wingen.name , ' window.'])
            t = (winlen/2+1):step:(length(obj.xx)-winlen/2-1);
            f = (0:dftn-1)/dftn*obj.fs;
            s = zeros(length(f) , length(t));
            for i=1:length(t)
                s(:,i) = fft(obj.xx(t(i)-winlen/2:t(i)+winlen/2-1) .* window , dftn);
            end
            t = t / obj.fs;
            
            if(isreal(obj.xx))
                % spectrum is symmetric -> only lower half is used and
                % presented
                f = f(1:ceil(length(f)/2));
                s = s(1:length(f),:);
            end
            
            if(nargin==6)
                figure(fignum)
%                 p1 = subplot(2,1,1);
                surf(t , f , 10*log10(s.*conj(s)) , 'EdgeColor','none');
                axis xy; axis tight; view(0,90);
                title(['Spectral Analysis using a ' , wingen.name , ' window']);
                xlabel('ms')
                ylabel('kHz')
%                 p2 = subplot(2,1,2);
%                 surf(t , f , angle(s)/pi , 'EdgeColor','none');
%                 axis xy; axis tight; view(0,90);
%                 linkaxes([p1,p2] , 'x');
            end
        end
    end
    
end

