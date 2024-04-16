% Class to control the Harvard Apparatus Pico Plus Pump. (Yet to add some advanced functions and the syringe database)
% #V.1-No errors
% All Volume Inputs should be in ul
% All Time Inputs in seconds
% All Lenths in mm

classdef PicoPlusElite
    properties(Access = private)
        BAUD = 115200; % Baud rate
        BITS = 8; % Data bits
        PARITY = 'none';
        STOP = 2; % Stop bits
        FC = 'none'; % Flow control
        Term = 'CR'; % Terminator
        minVelocity = 6e-5; % minimum pump velocity in [mm/s]
        maxVelocity = 1.18; % minimum pump velocity in [mm/s]
    end
    
    properties
        COM = [];
        verbose = 0;
    end
    
    properties(Dependent)
        diameter % Diameter of syringe in [mm]
        infuseRate % Infusion rate in [uL/s]
        withdrawRate % Withdraw rate in [uL/s]
        targetVolume %Target Volume in [ul]
        targetTime %Target time in [s]        
    end
    
    properties(Access = private)
        CR = char(13); % carriage return in ASCII
        SP = char(32); % space in ASCII
        LF = char(10); % line feed in ascii
    end

    methods % Basic methods
        
        function obj = PicoPlusElite(port)
            % This function creates an instance and establishes
            % communication with the hardware at the COM port specified.
            % If n=1, you get graphical assistance to choose the syringe.
            obj.COM = instrfind('Type','serial','Port',port,'Tag','');
            
            % Establish connection to instrument or clean up the existing
            % connection if the connection wasn't properly closed
            if isempty(obj.COM)
                obj.COM = serial(port);
            else
                fclose(obj.COM);
                obj.COM = obj.COM(1);
            end

            % Set serial communication settings
            set(obj.COM,'BaudRate',obj.BAUD,'DataBits',obj.BITS,...
               'Parity',obj.PARITY,'StopBits',obj.STOP,...
               'FlowControl',obj.FC,'Terminator',obj.Term)
            
            if strcmp(obj.COM.Status,'closed')
                fopen(obj.COM);
            end
        end
        
        function delete(obj)
            % Close the communication object and clean up
            fclose(obj.COM);
        end

        function Result = handler(obj)
            % Use this function to pass read back the status and response
            % from the pump after issuing a command
            pause(0.05)
            Result = fread(obj.COM,obj.COM.BytesAvailable);
            ind = find(Result == 13,1,'last');
            Status = char(Result(ind+2:end));
            Result = char(Result(2:ind-1)');
            if obj.verbose
                disp(Status)
                disp(Result)
            end
        end 
    end
    
    methods % Basic methods for communicating with hardware (ported functions)
        function Result = getAddress(obj)
            % Reads the address of the pump
            fprintf(obj.COM,'address');
            Result = obj.handler();
        end
        
        function Result = get.diameter(obj)
            % Displays the diameter set in the pump
            fprintf(obj.COM,'diameter');
            Result = obj.handler();
            Result = str2double(Result(1:end-2));
        end
        
        function obj = set.diameter(obj,d)
            %Sets the Diamter
            fprintf(obj.COM,['diameter ' num2str(d)]);
            obj.handler();
        end
          
        function Result = get.infuseRate(obj)
            % Reads the current Infusion Rate
            fprintf(obj.COM,'irate');
            Result = obj.handler();
            Result = str2double(Result(1:end-6));
        end
        
        function obj = set.infuseRate(obj,rate)
            % Sets the current infusion rate and makes sure it is a valid
            % value
            rate = obj.check_rate(rate);
            fprintf(obj.COM,['irate ' num2str(rate) ' ul/s']);
            obj.handler();
        end
        
        function Result = get.withdrawRate(obj)
            % Reads the current Withdraw Rate
            fprintf(obj.COM,'wrate');
            Result = obj.handler();
            Result = str2double(Result(1:end-6));
        end
        
        function obj = set.withdrawRate(obj,rate)
            % Sets the Withdraw Rate
            rate = obj.check_rate(rate);
            fprintf(obj.COM,['wrate ' num2str(rate) ' ul/s']);
            obj.handler();
        end
               
        function Result = get.targetVolume(obj)
            % Reads the current Target Volume
            fprintf(obj.COM,'tvolume');
            Result = obj.handler();
            Result = str2double(Result(1:end-2));
        end
        
        function obj = set.targetVolume(obj,v)
            % Sets the Target Volume
            fprintf(obj.COM,['tvolume ' num2str(v) ' ul']);
            obj.handler();
        end
        
        function Result = get.targetTime(obj)
            % Reads the current Target Time
            % NOTE: ramp time is also the target time and so this will also
            % modify the ramp time in the ramping mode
            fprintf(obj.COM,'ttime');
            Result = obj.handler();
            Result = str2double(Result(1:end-7));
        end
        
        function obj = set.targetTime(obj,v)
            % Sets the Target Time
            % NOTE: ramp time is also the target time and so this will also
            % modify the ramp time in the ramping mode
            fprintf(obj.COM,['ttime ' num2str(v)]);
            obj.handler();
        end
        
        function Result = get_infuseRamp(obj)
            % Gets the Infusion ramp settings:
            fprintf(obj.COM,'iramp');
            Result = obj.handler();
        end
        
        function obj = set_infuseRamp(obj,initialRate,finalRate,rampTime)
            % Sets the Infusion ramp profile.
            % initialRate = starting flow rate of ramp in [uL/s]
            % finalRate = ending flow rate of ramp in [uL/s]
            % rampTime = time for ramping from start to end rate in [s]
            % NOTE: ramp time is also the target time and this sets the
            % target time parameter simultaneously
            initialRate = obj.check_rate(initialRate);
            finalRate = obj.check_rate(finalRate);
            fprintf(obj.COM,['iramp ', num2str(initialRate), ' ul/s ',...
                 num2str(finalRate) ' ul/s ' num2str(rampTime)]);
            obj.handler();
        end 
        
        function Result = get_withdrawRamp(obj)
            %Gets the Withdraw ramp settings:
            fprintf(obj.COM,'wramp');
            Result = obj.handler();
        end

        function obj = set_withdrawRamp(obj,initialRate,finalRate,rampTime)
            % Sets the Withdraw ramp profile.
            % initialRate = starting flow rate of ramp in [uL/s]
            % finalRate = ending flow rate of ramp in [uL/s]
            % rampTime = time for ramping from start to end rate in [s]
            % NOTE: ramp time is also the target time and this sets the
            % target time parameter simultaneously
            initialRate = obj.check_rate(initialRate);
            finalRate = obj.check_rate(finalRate);
            fprintf(obj.COM,['wramp ', num2str(initialRate), ' ul/s ',...
                 num2str(finalRate) ' ul/s ' num2str(rampTime)]);
            obj.handler();
        end  

        function reset_volume(obj)
            %Resets the infused and withdrawn volumes 
            fprintf(obj.COM,'cvolume');
        end
        
        function reset_time(obj)
            %resets the infuse and withraw times
            fprintf(obj.COM,'ctime');
        end
        
        function run(obj)
        % Simulates a press of the run button
            fprintf(obj.COM,'run');
        end
        
        function stop(obj)
        % Instantly stops the pump
            fprintf(obj.COM,'stop');
        end
        
        function Result = get_running_rate(obj)
        % Displays the current rate that the motor is running at. A
        % valid response is returned only in dynamic situations
        % (while the pump is running). Quick Startmode only.
            fprintf(obj.COM,'crate');
            Result = obj.handler();
            Result = str2double(Result(1:end-6));
        end
        
        function Result = get_status(obj)
        % Gets the raw status for use with a controlling computer.
            fprintf(obj.COM,'status');
            Result = obj.handler();
        end
    end
    
    methods % custom helper functions
        function [minRate,maxRate] = get_rate_range(obj)
            % Use the current syringe diameter to calculate the minimum and
            % maximum flow rates that are possible.
            minRate = obj.minVelocity*pi/4*obj.diameter^2;
            maxRate = obj.maxVelocity*pi/4*obj.diameter^2;
        end
 
        function rate = check_rate(obj,rate)
            % Makes sure the rate is within the physically accesible range
            [minRate,maxRate] = obj.get_rate_range();
            if (rate < minRate) && (rate ~= 0);
                warning(['Flow rate out of range,',...
                    'setting to minimum rate: ', num2str(minRate), 'uL/s'])
                rate = minRate; 
            end
            if (rate > maxRate) && (rate ~= 0);
                warning(['Flow rate out of range,',...
                    'setting to maximum rate: ', num2str(maxRate), 'uL/s'])
                rate = maxRate; 
            end
        end
    end
    
    methods % Advanced methods for ease of use and custom functions        
        
        function set_rates(obj,rate)
            % Sets both the infuse rate and withdraw rate to the same value
            rate = obj.check_rate(rate);
            obj.withdrawRate = rate;
            obj.infuseRate = rate;
        end
        
        function infuse(obj, varargin)
        % Runs the pump in the infuse direction
        % Optional inputs are
        % rate in [uL/s]
        % targetVolume in [uL]
        obj.reset_volume();
            switch nargin
                case 1
                    fprintf(obj.COM,'irun');
                    obj.handler();
                case 2
                    obj.infuseRate = varargin{1};
                    fprintf(obj.COM,'irun'); 
                    obj.handler();
                case 3
                    obj.infuseRate = varargin{1};
                    obj.targetVolume = varargin{2};
                    fprintf(obj.COM,'irun'); 
                    obj.handler();
               case 4
                    obj.infuseRate = varargin{1};
                    obj.targetVolume = varargin{2};
                    fprintf(obj.COM,'irun'); 
                    obj.handler();
                    pause(obj.targetVolume/obj.infuseRate)
                otherwise
                    display('Invalid Arguments');
            end
        end
        
        function withdraw(obj, varargin)
        % Runs the pump in the infuse direction
        % Optional inputs are Rate followed by Target Volume
        obj.reset_volume();
            switch nargin
                case 1
                    fprintf(obj.COM,'wrun'); 
                case 2
                    obj.withdrawRate = varargin{1};
                    fprintf(obj.COM,'wrun'); 
                case 3
                    obj.withdrawRate = varargin{1};
                    obj.targetVolume = varargin{2};
                    fprintf(obj.COM,'wrun');
                case 4
                    obj.withdrawRate = varargin{1};
                    obj.targetVolume = varargin{2};
                    fprintf(obj.COM,'wrun'); 
                    obj.handler();
                    pause(obj.targetVolume/obj.withdrawRate)
                otherwise
                    display('Invalid Arguments');
            end
        end
        
        function infuse_and_withdraw(obj,rate,volume)
            obj.infuse(rate,volume,1);
            obj.withdraw(rate,volume,1);
        end
    end
end