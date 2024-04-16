clear all
baudrate = 4800;
port = 'COM9';
stopBits = 2;

s = serialport(port,baudrate, 'StopBits',stopBits);
configureTerminator(s,"CR")

char(read(s, s.NumBytesAvailable, "uint8"))
writeline(s,'?0')
%pause(0.01)
%R = char(read(s, s.NumBytesAvailable, "uint8"))
%pause(0.01)
%writeline(s,'?1')