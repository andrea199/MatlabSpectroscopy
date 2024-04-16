% Class to control the thor labs spectrophotometer 

% !!! Change instrument ID to the ID of your device !!!
%
% Type-ID:
% 0x8081 // CCS100 Compact Spectrometer
% 0x8083 // CCS125 Special Spectrometer 
% 0x8085 // CCS150 UV Spectrometer 
% 0x8087 // CCS175 NIR Spectrometer 
% 0x8089 // CCS200 UV-NIR Spectrometer 
%
% 'USB0::0x1313::<Type-ID>::<Serial Number>::0::RAW'
% Credits: http://materias.df.uba.ar/l5b2019c2/files/2013/02/Integrating-a-CCS-Spectrometer-in-MATLAB.pdf


classdef ThorSpectrophotometer
    properties(Access = private)
    libpath = 'C:\Program Files\IVI Foundation\VISA\Win64\';
    typeID = '0x8081';%'0x8089';  
    libname;  %TLCCS_64 file path
    hfile; %TLCCS.h file path
    statusA; % Device status codes 
    statusB; % Device status codes
    statusC; % Device status codes
    end
    
    properties 
        serialNo = 'M00803311';%'M00717497';
        inttime = 0.1; % Integration time (s)
        res; % Handles to the device
        hdl; % Handles to the device
        initError = 0;
    end
    
    methods % Basic functions
         function obj  = ThorSpectrophotometer(serialNo)
            %% Class constructor: Loads libraries and establishes connection
            if nargin==1
                obj.serialNo = serialNo;
            end
            obj.libname = strcat(obj.libpath, 'Bin\TLCCS_64.dll');
            obj.hfile = strcat(obj.libpath,'Include\TLCCS.h');
            
            loadlibrary(obj.libname,obj.hfile,'includepath',strcat(obj.libpath, 'Include\'), ...
                'includepath', strcat(obj.libpath,'Lib_x64\'));
            
            % Establish connection to the spectrophotometer
            obj.res=libpointer('int8Ptr',int8(strcat('USB0::0x1313::',obj.typeID,'::',obj.serialNo,'::0::RAW')));
            obj.hdl=libpointer('ulongPtr',0);
            [obj.statusA, obj.statusB, obj.statusC]=calllib('TLCCS_64', 'tlccs_init', obj.res, 0, 0, obj.hdl);
            
            % Display status if connectiion failed
            if obj.statusA~=0
                disp(['Initialize device failed with error code :', num2str(obj.statusA)]);
                obj.initError = 1;
            end   
            
        end   
        
        function delete(obj)
            % Class destructor
            calllib('TLCCS_64','tlccs_close', obj.hdl.value);
            unloadlibrary 'TLCCS_64';
            pause(1);
        end
        
        function seeLibraryCommands(obj)
            % See functions available in the loaded libarary
            libfunctionsview('TLCCS_64');
        end       
    end
    
    methods %Derived function
        function [wldata, specdata] = getSpectrum(obj)      
            
            %Get the complete spectrum
            calllib('TLCCS_64','tlccs_setIntegrationTime',obj.hdl.value,obj.inttime);
            calllib('TLCCS_64', 'tlccs_startScan', obj.hdl.value);            
            specdata=libpointer('doublePtr',double(1:3648));           
            calllib('TLCCS_64','tlccs_getScanData', obj.hdl.value, specdata);  
            wldata=libpointer('doublePtr',double(1:3648));
            calllib('TLCCS_64','tlccs_getWavelengthData', obj.hdl.value, 0, wldata, 0, 0);  
        end
        
        function binData  = getBinnedData(obj, binWavelength, binWidth)
            %Returns the binned data (average) for binWavelength, binWidth
            % binWavelength can be a nx1 vector, binWidth can be single
            % number of a a nx1 vector;
            % Note wlval and specval should be values from wldata and
            % specdata respectively
            
            numWavelenths = length(binWavelength);           
            if  numWavelenths~= length(binWidth)
                binWidth = binWidth*ones(1, numWavelenths);
            end
            
            %Acquire spectrum
            [wldata, specdata] = obj.getSpectrum();
            wlval = wldata.value; 
            specval = specdata.value;
                       
            % Calculate binData
            binData  = zeros(1, numWavelenths);
            for i = 1:numWavelenths
                [~, strtIdx] =  min(abs(wlval - (binWavelength(i) - binWidth(i))));
                [~, endIdx] = min(abs(wlval - (binWavelength(i) + binWidth(i))));
                binData(i) = mean(specval(strtIdx: endIdx));
            end    
        end
        function binData  = getBinnedDataHE(~,specval, strtIdx, endIdx)
            % High efficiency bin Data code. Assumes two wavelenths and
            % requires all the arguments
            binData  = [mean(specval(strtIdx(1): endIdx(1))), mean(specval(strtIdx(2): endIdx(2)))];
        end
        
    end
end