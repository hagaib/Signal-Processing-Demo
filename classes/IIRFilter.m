classdef IIRFilter < Filter
    %IIRFilter - Abstract class
    % Subclasses:
    % ButterFilter
    % Cheby1Filter
    % Cheby2Filter
    
    %NOTE: All filter implementations are not computationally efficient in
    %Matlab. They are implemented for educational purposes only
    
    properties
        isDigital % boolean (both analog and digital versions are supported)
        zeros
        poles
        gain
        Omegap  % passband angular frequencies (radians)
        Omegas  % stopband angular frequencies (radians)
    end
    
    methods (Access = public)
        function obj = IIRFilter(Fp , Fs , Ap , As , fs)
            
            obj@Filter(Fp , Fs);
            
            obj.Fp = Fp;
            obj.Fs = Fs;
            obj.Ap = Ap;
            obj.As = As;

            
            if(nargin < 5 || isempty(fs))
                obj.isDigital = false;
                %transform to angular frequencies:
                obj.Omegap = Fp*(2*pi);
                obj.Omegas = Fs*(2*pi);
            else
                obj.isDigital = true;
                %calculate normalized frequencies:
                obj.Wp = Fp/fs*2;
                obj.Ws = Fs/fs*2;
                %transform to angular frequencies using bilinear transform:
                obj.Omegap = tan(obj.Wp*pi/2);
                obj.Omegas = tan(obj.Ws*pi/2);
            end
            [z , p , g] = designLowpassPrototype(obj);
            [obj.zeros , obj.poles , obj.gain] = frequencyTransformation(obj , z , p , g);
            if(obj.isDigital)
                [obj.zeros , obj.poles , obj.gain] = IIRFilter.bilinearTransform(obj.zeros , obj.poles , obj.gain);
            end
            obj.N = max(length(obj.zeros) , length(obj.poles));
            [obj.b , obj.a] = IIRFilter.calculateAB(obj.zeros , obj.poles , obj.gain);
        end
        
            
        
        function plot(obj , fignum)
            try
                if(obj.isDigital)
                    sos = zp2sos(obj.zeros , obj.poles , obj.gain);
                    [h,w] = freqz(sos , 4096);
                else
                    error('go to catch');
                end
            catch
                if(obj.isDigital)
                    [h,w] = freqz(obj.b , obj.a , 4096);
                else
                    [h,w] = freqs(obj.b , obj.a , 4096);
                    w = w/2;
                end
            end
            figure(fignum)
            subplot(3,1,1)
            plot(w/(pi) , 10*log10(h.*conj(h)))
            title([class(obj) , ' magnitude response'])
            subplot(3,1,2)
            plot(w/(pi) , unwrap(angle(h)))
            title('phase response')
            subplot(3,1,3)
            plot(complex(obj.zeros) , 'o')            
            hold on
            plot(complex(obj.poles) , 'x' )
            if(obj.isDigital)
                plot(cos(0:0.1:2*pi) , sin(0:0.1:2*pi));
                plane = 'Z'; 
            else
                yL = ylim;
                line([0 0], yL); 
                plane = 'S'; 
            end
            hold off
            title(['zeros and poles on the ' , plane , ' plane.']);
            legend('zeros', 'poles');
        end
        
        
        %% Filter Implementation 2:
        %  Direct Form 2
        %%
        function sigout = filter2(obj , sigin)
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
            w = zeros(Ly,1);
            
            % Feed Back:
            for n = K+1:Ly
                w(n) = - obj.a(2:end) * w(n-1:-1:n-Na) + xx(n);
            end
            
            yy = zeros(Ly , 1);
            
            % Feed Forward:
            for n = K+1:Ly
                yy(n) =  obj.b * w(n:-1:n-Nb);
            end
            
            % Remove K zeros added at the beginning.
            yy = yy(K+1:Ly);
            
            
            if(isequal(class(sigin) , 'Signal'))
                sigout.xx = yy;
            else
                sigout = yy;
            end
        end
        
        %% Filter Implementation 3:
        %  Cascading Second Order Sections
        %%
        function sigout = filter3(obj , sigin)
            if(isequal(class(sigin) , 'Signal'))
                sigout = copy(sigin);
                sigout.title = [sigin.title , ' filtered by ' , class(obj)];
                xx = sigout.xx;
            else
                xx = sigin;
            end

            
            [sos , g] = zp2sos(obj.zeros , obj.poles , obj.gain);
            
            soscount = size(sos,1);
            yy = xx;
            for i=1:soscount
                bi = sos(i , 1:3);
                ai = sos(i , 4:6);
                cascadef = copy(obj);
                cascadef.a = ai;
                cascadef.b = bi;
                yy = cascadef.filter2(yy);
            end
            
            yy = g*yy;
            
            if(isequal(class(sigin) , 'Signal'))
                sigout.xx = yy;
            else
                sigout = yy;
            end
        end
        
        %% Filter Implementation 4:
        %  Parallel Form
        %%
        function sigout = filter4(obj , sigin)
            if(isequal(class(sigin) , 'Signal'))
                sigout = copy(sigin);
                sigout.title = [sigin.title , ' filtered by ' , class(obj)];
                xx = sigout.xx;
            else
                xx = sigin;
            end

            [r , p , C] = residuez(obj.b, obj.a);
            % FIR part:
            dummyfilt = copy(obj);
            dummyfilt.a = 1;
            dummyfilt.b = C;
            yy = dummyfilt.filter2(xx);
                        
            % IIR parallel filters:
            % Check if 1 order section is present
            if(mod(length(p), 2))
                [dummyfilt.b , dummyfilt.a] = residuez(r(end) , p(end) , []);
                yy = yy + dummyfilt.filter2(xx);
            end
            
            % 2nd order sections
            soscount = floor(length(p)/2);
            for i=1:soscount
                [bi , ai] = residuez(r((i*2-1):(i*2)) , p((i*2-1):(i*2)) , []);
                bi = real(bi);
                ai = real(ai);
                dummyfilt.a = ai; dummyfilt.b = bi;
                yy = yy + dummyfilt.filter2(xx);
            end
            
            if(isequal(class(sigin) , 'Signal'))
                sigout.xx = yy;
            else
                sigout = yy;
            end
        end
        
        
        %% Specifications of prototype analog low pass filter from the input specifications
        function [OmegaProP , OmegaProS] = prototypeLowPassSpecs(obj)
            if( isequal(obj.freqResType , 'lowpass') )
            % low pass to low pass
                OmegaProP = 1;
                OmegaProS = obj.Omegas/obj.Omegap;
            
            % low pass to high pass
            elseif( isequal(obj.freqResType , 'highpass') )
                OmegaProP = 1;
                OmegaProS = obj.Omegap/obj.Omegas;
            
            % low pass to band pass
            elseif( isequal(obj.freqResType , 'bandpass') )
                Bw = obj.Omegap(2) - obj.Omegap(1);
                Omega0_sq = obj.Omegap(1)*obj.Omegap(2);
                
                % Apply Geometric Symmetry:
                if( obj.Omegas(1) * obj.Omegas(2) < Omega0_sq )
                    % Make lower transition band tighter
                    obj.Omegas(1) = Omega0_sq / obj.Omegas(2);
                else
                    % Make higher transition band tighter
                    obj.Omegas(2) = Omega0_sq / obj.Omegas(1);
                end
                OmegaProP = 1;
                OmegaProS = (Omega0_sq - obj.Omegas(1)^2)/ (obj.Omegas(1) * Bw);
            end
        end
        
        
        %% Frequency transformations of zeros, poles and gain of analog prototype low pass filter according to frequency fresponse type
        function [z , p , g] = frequencyTransformation(obj , az , ap , ag)
            n = length(az);
            m = length(ap);
            if(isequal(obj.freqResType, 'lowpass'))
                z = az * obj.Omegap;
                p = ap * obj.Omegap;
                g = ag * obj.Omegap ^ (m-n);
            elseif(isequal(obj.freqResType, 'highpass'))
                z = obj.Omegap ./ az;
                p = obj.Omegap ./ ap;
                g = real(ag * (-1)^(n-m) * prod(az) / prod(ap));
                z = [z ; zeros(m-n , 1)];
                
            elseif(isequal(obj.freqResType, 'bandpass'))
                Bw = obj.Omegap(2) - obj.Omegap(1);
                Omega0_sq = prod(obj.Omegap);
                
                z = zeros(2*length(az) , 1);
                i=1;
                for root = az'
                    z(i:i+1) = roots([1 , -root*Bw , Omega0_sq]);
                    i=i+2;
                end
                z = [z ; zeros(m-n , 1)];
                
                p = zeros(2*length(ap) , 1);
                i=1;
                for root = ap'
                    p(i:i+1) = roots([1 , -root*Bw , Omega0_sq]);
                    i=i+2;
                end
                
                g = ag * Bw ^(m-n);
            end
        end
        
    end
    
    
    methods (Access = protected)
        function [z,p,g] = designLowpassPrototype(obj)
            % To be implemented by subclasses
        end

    end
    
        
    methods (Static)
        
        function [z , p , g] = bilinearTransform(zeros , poles , gain)
            % transfers zeros poles and gain of analog filter to discrete
            % filter using bilinear transform
            p = (1+poles)./(1-poles);
            z = (1+zeros)./(1-zeros);
            addZeros = length(poles) - length(zeros);
            z = [z ; -ones(addZeros , 1)];
            g = real(gain * prod(1- zeros) / prod(1-poles));
        end
        
        
        function [b , a] = calculateAB(zeros , poles , gain)
            b = real(gain * poly(zeros));
            a = real(poly(poles));
        end
        
    end
    
end

