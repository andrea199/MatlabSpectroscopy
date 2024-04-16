% Class to control the Harvard Apparatus Pico Plus Pump. (Yet to add some advanced functions and the syringe database)
% #V.1-No errors
% All Volume Inputs should be in ul
% All Time Inputs in seconds
% All Lenths in mm

classdef Ainsworth
    properties(Access = private)
        BAUD = 4800; % Baud rate
        BITS = 8; % Data bits
        PARITY = 'none';
        STOP = 2; % Stop bits
        FC = 'none'; % Flow control
        Term = 'CR'; % Terminator
    end
    
    properties
        COM = [];
        verbose = 0;
    end
    properties(Dependent)
        mass % Mass in [g]    
    end
     properties(Access = private)
        CR = char(13); % carriage return in ASCII
        SP = char(32); % space in ASCII
        LF = char(10); % line feed in ascii
     end
     methods % Basic methods
        
        function obj = Ainsworth(port)
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
               'FlowControl',obj.FC,'Terminator',obj.Term, 'InputBufferSize', 2048*2048)
            
            if strcmp(obj.COM.Status,'closed')
                fopen(obj.COM);
                fprintf(obj.COM,'?0');
            end
            pause(0.1)
        end
        
        function delete(obj)
            % Close the communication object and clean up
            fclose(obj.COM);
        end

        function Result = handler(obj)
            % Use this function to pass read back the status and response
            % from the pump after issuing a command

             Result = fread(obj.COM,obj.COM.BytesAvailable);
             ind = find(Result == 43,1,'first');
             ind2 = find(Result(ind:end) == 13,1,'first');

           %  Status = char(Result(ind+2:end));
             Result = char(Result(ind+1:ind+ind2-1)');
            if obj.verbose
              %  disp(Status)
                disp(Result)
             end
        end 
     end
   methods % Basic methods for communicating with hardware (ported functions)
        function Result = get.mass(obj)
            % Displays the mass immediately
            if obj.COM.BytesAvailable>30
                fread(obj.COM,obj.COM.BytesAvailable-30);
            end
            try            
            Result = obj.handler();
            Result = str2double(Result(1:end));
            catch
                Result=nan;
            end
        end 
   end
end
