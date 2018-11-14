classdef Component < handle
    
    properties
        I_SIGNAL_CODE   % 1=white noise , 2=CW (sinusoid) 3=PMCW (pulse modulated CW)
        A               % amplitude (or stdev for noise component)
        F               % frequency (khz)
        TAU             % pulse duration
        PRI             % pulse repetition interval
        N               % number of pulses 
        I_PHASE_CODE    % 1=every pulse starts with a different phase (all pulses are generated using the same oscillator). 2=every pulse starts with the same phase
        xx              % signal samples
        fs              % sample rate (k samples/sec)
    end
    
    methods
        function obj = Component(arr)
            if(nargin<1)
                error('Wrong number of parameters');
            end
            pcount = length(arr);
            if(pcount ~= 7)
                if (arr(1)==1 && pcount ~= 2)
                    error(['Component of code ' , num2str(arr(1)) , ' must have 2 parameters']);
                elseif(arr(1)==2 && pcount ~= 3)
                    error(['Component of code ' , num2str(arr(1)) , ' must have 3 parameters']);
                elseif(arr(1)==3 && pcount ~= 7)
                    error(['Component of code ' , num2str(arr(1)) , ' must have 7 parameters']);
                end
            end
            if (arr(1) < 1 || arr(1) > 3)
                error('Wrong signal code. Possible values are: 1,2,3');
            end
            obj.I_SIGNAL_CODE = arr(1);
            obj.A = arr(2);
            if(obj.I_SIGNAL_CODE ~= 1)
                obj.F = arr(3);
            end
            if(obj.I_SIGNAL_CODE == 3)
                obj.TAU = arr(4);
                obj.PRI = arr(5);
                obj.N = arr(6);
                obj.I_PHASE_CODE = arr(7);
            end
        end
        
        function xx = generate(obj, sample_rate , time)
            %time in seconds
            %sample rate in hz
            obj.fs = sample_rate;
            switch (obj.I_SIGNAL_CODE)
                case 1
                    obj.xx = gaussian_sample(obj , floor(sample_rate*time));
                case 2
                    obj.xx = CW_sample(obj , sample_rate , time);
                case 3
                    obj.xx = pulse_sample(obj , sample_rate,time);
                otherwise
                    
            end
            xx=obj.xx;
        end
        
    end
    
    
    
    methods(Access = protected)
        function x = pulse_sample(obj , sample_rate , time)
            if(obj.I_PHASE_CODE==1)
                 x = CW_sample(obj , sample_rate , obj.PRI * obj.N);
                 for k=0:obj.N-1
                     x(floor((k*obj.PRI+obj.TAU)*sample_rate)+1:floor((k+1)*obj.PRI*sample_rate))=0;
                 end
                 
            elseif(obj.I_PHASE_CODE==2)
                T_samples = floor(obj.PRI * sample_rate);
                one_rep = zeros(T_samples , 1);
                one_rep(1:floor(obj.TAU*sample_rate)) = ...
                    CW_sample(obj , sample_rate, obj.TAU);
                x = zeros(T_samples * obj.N , 1);
                for k=0:obj.N-1
                     x(k*T_samples+1:(k+1)*T_samples) = one_rep;
                 end
            end
            
            sample_count = floor(sample_rate*time);
            diff = sample_count - length(x);
            if(diff >= 0)
%                 x = [x ; zeros(floor(diff),1)];
                x = [zeros(floor(diff/2),1); x; zeros(ceil(diff/2) , 1)];
            else
                error('time should be at least: %0.2f seconds\n' , length(x)/sample_rate);
            end
            
        end

        function x = CW_sample(obj, sample_rate , time)
            k = floor(sample_rate * time);
            cycle = sample_rate / obj.F;
            x = obj.A * cos(2*pi/cycle*(0:(k-1)));
            x = transpose(x);
        end

        function x = gaussian_sample(obj,k)
            %using box-muller algorithm
            x = sqrt(-2.0*log(rand(k,1))).*cos(2*pi*rand(k,1));
            x = obj.A * x;
        end

    end
end