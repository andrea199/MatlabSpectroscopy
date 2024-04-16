classdef ArduinoPump
    properties(Access = private)
        BAUD = 9600; % Baud rate
        BITS = 8; % Data bits
        PARITY = 'none';
        STOP = 2; % Stop bits
        FC = 'none'; % Flow control
        Term = 'LF'; % Terminator
    end
    
    properties
        COM = [];
        verbose = 0;
    end
    properties(Dependent)
    end
     properties(Access = private)
        CR = char(13); % carriage return in ASCII
        SP = char(32); % space in ASCII
        LF = char(10); % line feed in ascii
     end
     methods % Basic methods
        
        function obj = ArduinoPump(port)
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
        
        function actuateMotor(obj,str)
            % Reads the address of the pump
            fprintf(obj.COM,str);
        end
     end
end