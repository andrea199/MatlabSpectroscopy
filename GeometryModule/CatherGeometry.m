classdef CatherGeometry < handle
    % Class for retreiving the geometry features of catheters necessary for computations based
    % on the make and model, and the injection ports used for drug
    % dleivery. 
    % P1 - Closest to patient
    properties
        make;  %Make of the catheter
        model; %Model of the cather
        ports; %Ports with p1 being the port closest to the patient. Lenth of vector will be equal to numDrugs  
        params; %Struture to hold the geometry based on the make and model
        nGrid; % delta x for spatial discretization
        geometry; % Structure with cell with radius vectors, domain length. 
        radSmoothVal = 1; %Radii smooth param
    end
    
    methods % Constructor
        function obj = CatherGeometry(make, model)
          obj.make = make;
          obj.model = model;         
        end           
    end 
    
    methods % Methods to calculate the radii ofthe catheters
        
        function parseGeometry(obj,ports, nGrid)
            %v1: Only supports single nGrid (no adaptive or spatially variable nGrid)
          obj.ports = sort(ports); 
          obj.nGrid = nGrid;
          switch strcat(obj.make, obj.model)
              case 'BBN456003' %All units in mm
                  obj.params.p1toTubing = 30; %Port 1 to start of tubing % Changed from 22 based on CT scan
                  obj.params.ptop  = 34;   % Distance between ports % Changed from 33
                  obj.params.tubing2Junction = 104; % Tubing to Junction multilumen junction
                  obj.params.junction2Patient = 190; % Multilumen junction to catheter tip
                  obj.params.radii.stopcock = 2; %Inner average radii of the stopcock %Changed to from 2.1
                  obj.params.radii.tubing = 0.75; %Inner average radii of the tubing % Changed from 0.8 based on CT scan
                  obj.params.radii.multilumenCatheter =0.5; %Inner average of the multilumen catheter  %Changed to 0.5 based on CT scan image
                  obj.params.ghostDomainLength = 50; % mm 
                  obj.calculateRadii()
              case 'Cylinder' %All units in mm, Test case for Gill-Subramian validation
                  obj.params.radii.tubing = 1;                  
                  
                  %Perform calculations equivalent to calculateRadii()
                  % 11/05/22 - Added realDomainLength & realDomainIndex to
                  % account for a ghost section to avoid finite boundary effects
                  % qcorr will be truncated by the realDomainIndex
                  
                  obj.geometry.realDomainLength = obj.ports;
                  obj.geometry.realDomainIndex = obj.ports/obj.nGrid + 1;
                  obj.geometry.domainLength(1) = obj.ports + 25; %Ghost domain length of 25                  
                  obj.geometry.radii{1} = ones(length(0:obj.nGrid:obj.geometry.domainLength(1)),1)*obj.params.radii.tubing/1000;
             
              case 'StepCylinder' %All units in mm, Ports specify length of the ports with radii of 2mm for the first 
                  %half and 1 mm for the second half. 
                  obj.params.radii.tubing = 1;                  
                  
                  %Perform calculations equivalent to calculateRadii()
                  % 11/05/22 - Added realDomainLength & realDomainIndex to
                  % account for a ghost section to avoid finite boundary effects
                  % qcorr will be truncated by the realDomainIndex
                  
                  obj.geometry.realDomainLength = obj.ports;
                  obj.geometry.realDomainIndex = obj.ports/obj.nGrid + 1;
                  obj.params.ghostDomainLength = 25;
                  obj.geometry.domainLength(1) = obj.ports + obj.params.ghostDomainLength; %Ghost domain length of 25                  
                  obj.geometry.radii{1} = smooth([ones(1,length(0:obj.nGrid:(obj.geometry.realDomainLength(1)/2)-obj.nGrid))*2*obj.params.radii.tubing/1000, ...
                                          ones(1, length(0:obj.nGrid:obj.geometry.realDomainLength(1)/2+obj.params.ghostDomainLength))*obj.params.radii.tubing/1000],obj.radSmoothVal);   
                  
                  
          end          
          
        end  
        
        function calculateRadii(obj)
            % Calculates the radii based on params, ports and nGrid
            % settings                
                obj.geometry.firstDrug2Tubing = (obj.ports(1)-1)*obj.params.ptop + obj.params.p1toTubing; % Distance of first drug to tubing
                obj.geometry.domainLength(1) = obj.geometry.firstDrug2Tubing + obj.params.tubing2Junction + (obj.params.junction2Patient + obj.params.ghostDomainLength); %Ghost domain length of 50 units with radii 
                obj.geometry.radii{1} = smooth([ones(1, length(0:obj.nGrid:obj.geometry.firstDrug2Tubing - obj.nGrid))*obj.params.radii.stopcock/1000,...
                                         ones(1,length(0:obj.nGrid:obj.params.tubing2Junction-obj.nGrid))*obj.params.radii.tubing/1000,...
                                         ones(1,length(0:obj.nGrid:(obj.params.junction2Patient + obj.params.ghostDomainLength)))*obj.params.radii.multilumenCatheter/1000;
                                        ],obj.radSmoothVal);
                                     
                obj.geometry.realDomainLength(1) = obj.geometry.domainLength(1) - obj.params.ghostDomainLength;
                obj.geometry.realDomainIndex = length(obj.geometry.radii{1}) - obj.params.ghostDomainLength/obj.nGrid;
            
            % Pending: Domain length correction for multi-drug systems
            
            for i = 2:length(obj.ports)
                obj.geometry.domainLength(i) = (obj.ports(i)-obj.ports(i-1))*obj.params.ptop;
                obj.geometry.radii{i} = ones(length(0:obj.nGrid:obj.geometry.domainLength(i)-obj.nGrid),1)*obj.params.radii.stopcock/1000;
            end        
        end
    
    end
end
