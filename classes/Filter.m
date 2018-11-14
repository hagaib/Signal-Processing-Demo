classdef Filter < matlab.mixin.Copyable
    % Abstract 
    % Superclass for all frequency selective filter classes
    % Only lowpass , highpass and bandpass filters are supported
    % NOTE: all filtering implementations are not efficient in Matlab.
    % They are implemented for educational purposes only.
    
    
    properties
        freqResType % values: 'lowpass' , 'highpass' , 'bandpass'.
        N           % filter order
        Fp          % passband frequency (hertz)
        Fs          % stopband frequency (hertz)
        Ap          % passband amplitude (db)
        As          % stopband amplitude (db)
        Wp          % passband normalized frequency
        Ws          % stopband normalized frequency
        b           % transfer function feedforward coefficients
        a           % transfer function feedback coefficients
    end
    
    methods
        % abstract constructor
        % lowpass or highpass filters: Fp , Fs are scalar.
        % bandpass filters: Fp , Fs are 2 dimensional vectors
        function obj = Filter(Fp , Fs)
            if( length(Fp) >1 && length(Fs) >1)
                if ( Fp(1)>Fs(1) && Fp(2) < Fs(2))
                    obj.freqResType = 'bandpass';
                else
                    error('bandstop not supported.');
                end
            elseif(length(Fp) == 1 && length(Fs) == 1)
                if(Fp < Fs)
                    obj.freqResType = 'lowpass';
                elseif(Fp > Fs)
                    obj.freqResType = 'highpass';
                end
                    
            else
                error('Frequency response specified is not supported.');
            end
            
            disp([obj.freqResType , ' filter of class ' class(obj) , ' is generated.']);
        end
        
        function plot(obj , fs , figurenum)
            if(nargin < 2)  , fs = []; end
            if(nargin <3)
                figure
            else
                figure(figurenum)
            end
            
            if(isempty(fs))
                freqz(obj.b , obj.a , 1000);
            else
                freqz(obj.b , obj.a , 1000 , fs);
            end
        end
        
        
        % Filter Implementations: Direct Form 1
        % Other implementations are irrelevant to FIR Filtering, so their
        % implementation is present in IIRFilter class.
        function sigout = filter1(obj , sigin)
            if(isequal(class(sigin) , 'Signal'))
                sigout = copy(sigin);
                sigout.title = [sigin.title , ' filtered by ' , class(obj)];
                xx = sigout.xx;
            else
                xx = sigin;
            end

            Nb = length(obj.b)-1; 
            Na = length(obj.a)-1; 
            K = max(Na,Nb); 
            Lx = length(xx);
            % add K samples to beginning of signal so no need to worry
            % about transients
            xx = [zeros(K,1) ; xx(:)];
            Ly = Lx+K; 
            yy = zeros(Ly,1);
            
            % Feed Forward:
            for n = K+1:Ly
                yy(n) = obj.b*xx(n:-1:n-Nb);
            end
            
            % Feed Back:
            for n = K+1:Ly
                yy(n) = yy(n) - obj.a(2:end) * yy(n-1:-1:n-Na);
            end
            
            % Remove K zeros added at the beginning.
            yy = yy(K+1:Ly);
            
            
            if(isequal(class(sigin) , 'Signal'))
                sigout.xx = yy;
            else
                sigout = yy;
            end
        end
    end
    
end

