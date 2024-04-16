classdef MultiDrugSolver < handle
    %Class with routines to solve multidrug profiles based on Solvers
    %defined in ADEModule
    
    properties
        drugDiffusivity; %In v1 all drugs have same diffusivity
        DeltaT = 1; %In v1 the solver has fixed time step and cannot be externally controlled
        ADE; %Instances of the ADEModule
        ICM; %Instance of the IntegralControlModule
        CG; %Instance of the CatheterGeometry
        Drug = {}; %Holds the data on the drug concentration
        % Post Python updates
        lamdN; % Holds the zeros of J1 (bessel function of the first kind)
        B1n; % Coefficients of Gill-Sub Diffusion coefficient
        currTime = 1; % Current cumulative time
        %Scale = []; % Scale factor for correcting the concentration
        %PumpFlowRatesPrevStep; % Pump flows in the previous time step.t
    end
    
    methods % Constructor
        function obj = MultiDrugSolver(CG, ICM)
          %Initialize ADE Modules based on data from the CGModule     
          for i = 1:ICM.Drugs.NumDrugs
              obj.ADE{i} = ADEModule(); 
              obj.ADE{i}.Ngrid = length(CG.geometry.radii{i});
              obj.ADE{i}.RoLScale = CG.geometry.radii{i};
              obj.ADE{i}.DomainLength = CG.geometry.domainLength(i)/1000; %Convert to meters    
              obj.ADE{i}.RealDomainIndex = CG.geometry.realDomainIndex;
              obj.ADE{i}.DcoeffFn = @(u) obj.drugDiffusivity*(1 + ((u.*obj.ADE{i}.RoLScale/obj.drugDiffusivity).^2)/48);
          end
          %Gill-Sub Diffusion coefficient calculation
          obj.lamdN = obj.besselzero(1, 10, 1);
          obj.B1n = besselj(3,obj.lamdN).* besselj(2,obj.lamdN)./(( besselj(0,obj.lamdN).^2).*obj.lamdN.^5);
          
          obj.ICM = ICM;
          obj.CG = CG;
        end    
    end
    
    methods % Core functions
        
        function obj = CalculateNextTimeStep(obj)
            %Evolves the multi drug concentration by a single time step  
            obj.initializeSolvers();
            
            %Solve the primary drug concentration (at the injection area)
            % Solve from upstream
            % Scale appropriately (v2), Generalize formultiple drugs
            for i = 1:obj.ICM.Drugs.NumDrugs 
                obj.ADE{i}.Solver.q = obj.Drug{i,i}; %Update the concentration from the previous time step (Scaling needs to be applied here)
                obj.ADE{i}.Solver.q(1) = (obj.ICM.Pumps.PumpFlowRate(i+1)/(obj.ICM.Pumps.PumpFlowRate(1)+obj.ICM.Pumps.PumpFlowRate(i+1))); %**Fix error** in the calculation of upstream drug flux. 
                obj.ADE{i}.stepCrankNicolson();     
                obj.Drug{i,i} = obj.ADE{i}.Solver.q;
                %obj.ADE{i}.Solver.qcorr(:,i) = obj.Drug{i,i}(1:obj.ADE{i}.RealDomainIndex,1).*((obj.ADE{i}.RoLScale(1)^2)./obj.ADE{i}.RoLScale(1:obj.ADE{i}.RealDomainIndex).^2)/obj.Drug{i,i}(1,1); %Removed normalization based on new definition            
                obj.ADE{i}.Solver.qcorr(:,i) = obj.ADE{i}.Solver.q/obj.ADE{i}.Solver.q(1); % After correction of M and N, normalization with RoLScale is not required
            end     
           
           %Pending multidrug changes: (a) Domain truncation 
           % Down stream propogation 
           for i = 2:obj.ICM.Drugs.NumDrugs
                for j=i-1:-1:1 
                obj.Drug{i,j}(1,1) = obj.Drug{i,j+1}(end,1); %Matching upstream drug 
                obj.ADE{j}.Solver.q = obj.Drug{i,j};   
                obj.ADE{j}.stepCrankNicolson(); 
                obj.Drug{i,j} = obj.ADE{j}.Solver.q ;  
                %obj.ADE{j}.Solver.qcorr(:,i) = obj.Drug{i,j}.*((obj.ADE{j}.RoLScale(1)^2)./obj.ADE{j}.RoLScale.^2)/obj.Drug{i,i}(1,1); %Removed normalization based on new definition            
                obj.ADE{j}.Solver.qcorr(:,i) = obj.ADE{i}.Solver.q/obj.ADE{i}.Solver.q(1);
                end               
           end 
            
            
        end       
        
        function obj = CalculateFirstTimeStep(obj)
            %Calculates the first time step (strong gradients). Currently, the code only
            % works for t = 0. However with modification the same function can be used when 
            %a new drug is impulsively introduced into the system at an
            %arbitary time
            %              Drug 2, P2            Drug 1, P1       
            %-----------------|  |-------------------|  |-----------
            %-------------------------------------------------------
            
            % Initialize drug concentrations
            for i = 1:obj.ICM.Drugs.NumDrugs
                obj.Drug(i:obj.ICM.Drugs.NumDrugs,i) = mat2cell(zeros(length(obj.ADE{i}.RoLScale),obj.ICM.Drugs.NumDrugs-i+1),length(obj.ADE{i}.RoLScale),ones(1,obj.ICM.Drugs.NumDrugs-i+1));
            end 
            
            % Initialize the CN solver
            obj.initializeSolvers();
            
            %Solve the primary drug concentration (via BE)
            for i = 1:obj.ICM.Drugs.NumDrugs              
                obj.ADE{i}.initializeBackwardEuler();                
                %Make concentration dimensional
                obj.ADE{i}.Solver.q = obj.ADE{i}.Solver.q*(obj.ICM.Pumps.PumpFlowRate(i+1)/(obj.ICM.Pumps.PumpFlowRate(1)+obj.ICM.Pumps.PumpFlowRate(i+1)));
                obj.ADE{i}.stepBackwardEuler();  
                obj.Drug{i,i} = obj.ADE{i}.Solver.q;
                %obj.ADE{i}.Solver.qcorr(:,i) = obj.Drug{i,i}(1:obj.ADE{i}.RealDomainIndex,1).*((obj.ADE{i}.RoLScale(1)^2)./obj.ADE{i}.RoLScale(1:obj.ADE{i}.RealDomainIndex).^2)/obj.Drug{i,i}(1,1); %Removed normalization based on new definition            
                obj.ADE{i}.Solver.qcorr(:,i) = obj.ADE{i}.Solver.q/obj.ADE{i}.Solver.q(1);
            end  
            
            %Solve drug concentration downstream (via BE). i - drug. j - position
            % Current state - Needs validation
            for i = 2:obj.ICM.Drugs.NumDrugs
                for j=i-1:-1:1 
                obj.Drug{i,j}(1,1) = obj.Drug{i,j+1}(end,1); %Matching upstream drug 
                obj.ADE{j}.Solver.q = obj.Drug{i,j};   
                obj.ADE{j}.stepBackwardEuler(); 
                obj.Drug{i,j} = obj.ADE{j}.Solver.q ;  
                %obj.ADE{j}.Solver.qcorr(:,i) = obj.Drug{i,j}.*((obj.ADE{j}.RoLScale(1)^2)./obj.ADE{j}.RoLScale.^2)/obj.Drug{i,i}(1,1);                 
                obj.ADE{j}.Solver.qcorr(:,i) = obj.ADE{j}.Solver.q/obj.ADE{j}.Solver.q(1);
                end               
            end
            
        end
    end
        
        
        methods (Access=private) %Supporting functions
            
            function obj = initializeSolvers(obj)
                % Initalize the CN solver on all modules (for all Drugs)             
                CarrierFlowRate = obj.ICM.Pumps.PumpFlowRate(1);
                for i = 1:obj.ICM.Drugs.NumDrugs
                    obj.ADE{i}.UScale = (CarrierFlowRate + sum(obj.ICM.Pumps.PumpFlowRate(i+1:end)))./(pi*obj.ADE{i}.RoLScale.^2); 
                    obj.ADE{i}.Dcoeff = feval(obj.ADE{i}.DcoeffFn, obj.ADE{i}.UScale);
                    %obj.DcoeffGillSub(i); % GillSub Diffcusion coefficient
                    obj.ADE{i}.DeltaT = obj.DeltaT;
                    obj.ADE{i}.initializeCrankNicolson();
                end
            end
            
            function obj = DcoeffGillSub(obj, ind)
                % Calculate the diffusion coefficient based on Gill-Sub
                % (1970). Pe is defined on the max velocity
                
                Pe = 2*obj.ADE{ind}.UScale.*obj.ADE{ind}.RoLScale/obj.drugDiffusivity;
                %tau = 0.5*obj.ADE{ind}.DomainLength./(obj.ADE{ind}.RoLScale.*Pe); % Remove the extra factor of 0.5
                tau = obj.drugDiffusivity*obj.currTime./(obj.ADE{ind}.RoLScale.^2);
                deltaD = sum(obj.B1n.*(exp(-(obj.lamdN.^2).*tau)),2);
                obj.ADE{ind}.Dcoeff = obj.drugDiffusivity*(1 + (Pe.^2)/192 -(4*Pe.^2).*deltaD);                
            end
            
            % Supporting functions
            function x = besselzero(obj, n, k, kind)
            % besselzero calculates the zeros of Bessel function of the first and second kind 
            %
            %   x = besselzero(n)
            %   x = besselzero(n, k)
            %   x = besselzero(n, k, kind)
            %
            %% Inputs
            % * *n* - The order of the bessel function. n can be a scalar, vector, or
            %       matrix.  n can be positive, negative, fractional, or any
            %       combinaiton. abs(n) must be less than or equal to
            %       146222.16674537213 or 370030.762407380 for first and second kind
            %       respectively. Above these values, this algorithm will not find the
            %       correct zeros because of the starting values therefore an error is
            %       thrown instead.
            % * k - The number of postive zeros to calculate.  When k is not supplied,
            %       k = 5 is the default. k must be a scalar.
            % * kind - kind is either 1 or 2. When kind is not supplied, default is
            %          kind = 1.
            %
            %% Outputs
            % * x - The calculated zeros.  size(x) = [size(n) k].
            % 
            %% Description
            % besselzero calculates the first k positive zeros of nth order bessel
            % function of the first or second kind.  Note, that zero is not included as
            % the first zero.
            %
            %% Algorithm
            % the first three roots of any order bessel can be approximated by a simple
            % equations.  These equations were generated using a least squares fit of
            % the roots from orders of n=0:10000. The approximation is used to start
            % the iteration of Halley's method.  The 4th and higher roots can be
            % approximated by understanding the roots are regularly spaced for a given
            % order.  Once the 2nd and 3rd roots are found, the spacing can be
            % approximated by the distance between the 2nd and 3rd root.  Then again
            % Halley's method can be applied to precisely locate the root.
            %%
            % Because the algorithm depends on good guesses of the first three zeros,
            % if the guess is to far away then Halley's method will converge to the
            % wrong zero which will subsequently cause any other zero to be incorrectly
            % located. Therefore, a limit is put on abs(n) of 146222.16674537213 and
            % 370030.762407380 for first and second kind respectively.  If n is
            % specified above these limits, then an error is thrown.
            %
            %% Example
            %   n = (1:2)';
            %   k = 10;
            %   kind = 1;
            %   z = besselzero(n, k, kind);
            %   x = linspace(0, z(end), 1000);
            %   y = nan(2, length(x));
            %   y(1,:) = besselj(n(1), x);
            %   y(2,:) = besselj(n(2), x);
            %   nz = nan(size(z));
            %   nz(1,:) = besselj(n(1), z(1,:));
            %   nz(2,:) = besselj(n(2), z(2,:));
            %   plot(x, y, z, nz,'kx')
            %
            % Originally written by 
            % Written by: Greg von Winckel - 01/25/05
            % Contact: gregvw(at)chtm(dot)unm(dot)edu
            %
            % Modified, Improved, and Documented by 
            % Jason Nicholson 2014-Nov-06
            % Contact: jashale@yahoo.com
            %% Change Log
            % * Original release. 2005-Jan-25, Greg von Winckel.
            % * Updated Documentation and commented algorithm. Fixed bug in finding the
            %   the first zero of the bessel function of the second kind. Improved speed by 
            %   factor of 20. 2014-Nov-06, Jason Nicholson.
            %
            %% Input checking
            assert(nargin >=1 | nargin <=3,'Wrong number of input arguments.');
            % Take care of default cases of k and kind
            if nargin < 2
                k = 5;
            end
            if nargin < 3
                kind = 1;
            end
            assert(isscalar(kind) & any(kind == [1 2]), '''kind''must be a scalar with value 1 or 2 only.');
            assert(isscalar(k) & fix(k)==k & k>0, 'k must a positive scalar integer.');
            assert(all(isreal(n(:))), 'n must be a real number.');
            % negative orders have the same roots as the positive orders
            n = abs(n);
            % Check for that n is less than the ORDER_MAX
            if kind==1
                ORDER_MAX = 146222.16674537213;
                assert(all(n <= ORDER_MAX), 'all n values must be less than or equal %6.10f for kind=1.', ORDER_MAX);
            elseif kind==2
                ORDER_MAX = 370030.762407380;
                assert(all(n(:) <= ORDER_MAX), 'all n values must be less than or equal %6.10f for kind=2.', ORDER_MAX);
            end
            %% Setup Arrays
            % output size
            nSize = size(n);
            if nSize(end) ==1
                outputSize = [nSize(1:end-1) k];
            else
                outputSize = [nSize k];
            end
            % number of orders for each kth root
            nOrdersPerRoot = prod(outputSize(1:end-1));
            x = nan(outputSize);
            %% Solve for Roots
            switch kind
                case 1
                    % coefficients and exponent are from least squares fitting the k=1,
                    % n=0:10000.
                    coefficients1j = [0.411557013144507;0.999986723293410;0.698028985524484;1.06977507291468];
                    exponent1j = [0.335300369843979,0.339671493811664];
                    % guess for k = 1
                    x((1:nOrdersPerRoot)') = coefficients1j(1) + coefficients1j(2)*n(:) + coefficients1j(3)*(n(:)+1).^(exponent1j(1)) + coefficients1j(4)*(n(:)+1).^(exponent1j(2));
                    % find first zero
                    x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));

                    if k >= 2
                        % coefficients and exponent are from least squares fitting the k=2,
                        % n=0:10000.
                        coefficients2j = [1.93395115137444;1.00007656297072;-0.805720018377132;3.38764629174694];
                        exponent2j = [0.456215294517928,0.388380341189200];
                        % guess for k = 2
                        x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2j(1) + coefficients2j(2)*n(:) + coefficients2j(3)*(n(:)+1).^(exponent2j(1)) + coefficients2j(4)*(n(:)+1).^(exponent2j(2));
                        % find second zero
                        x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
                    end

                    if k >= 3
                        % coefficients and exponent are from least squares fitting the k=3,
                        % n=0:10000.
                        coefficients3j = [5.40770803992613;1.00093850589418;2.66926179799040;-0.174925559314932];
                        exponent3j = [0.429702214054531,0.633480051735955];
                        % guess for k = 3
                        x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3j(1) + coefficients3j(2)*n(:) + coefficients3j(3)*(n(:)+1).^(exponent3j(1)) + coefficients3j(4)*(n(:)+1).^(exponent3j(2));
                        % find second zero
                        x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
                    end
                case 2
                    % coefficients and exponent are from least squares fitting the k=1,
                    % n=0:10000.
                    coefficients1y = [0.0795046982450635;0.999998378297752;0.890380645613825;0.0270604048106402];
                    exponent1y = [0.335377217953294,0.308720059086699];
                    % guess for k = 1
                    x((1:nOrdersPerRoot)') = coefficients1y(1) + coefficients1y(2)*n(:) + coefficients1y(3)*(n(:)+1).^(exponent1y(1)) + coefficients1y(4)*(n(:)+1).^(exponent1y(2));
                    % find first zero
                    x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));

                    if k >= 2
                        % coefficients and exponent are from least squares fitting the k=2,
                        % n=0:10000.
                        coefficients2y = [1.04502538172394;1.00002054874161;-0.437921325402985;2.70113114990400];
                        exponent2y = [0.434823025111322,0.366245194174671];
                        % guess for k = 2
                        x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2y(1) + coefficients2y(2)*n(:) + coefficients2y(3)*(n(:)+1).^(exponent2y(1)) + coefficients2y(4)*(n(:)+1).^(exponent2y(2));
                        % find second zero
                        x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
                    end

                    if k >= 3
                        % coefficients and exponent are from least squares fitting the k=3,
                        % n=0:10000.
                        coefficients3y = [3.72777931751914;1.00035294977757;2.68566718444899;-0.112980454967090];
                        exponent3y = [0.398247585896959,0.604770035236606];
                        % guess for k = 3
                        x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3y(1) + coefficients3y(2)*n(:) + coefficients3y(3)*(n(:)+1).^(exponent3y(1)) + coefficients3y(4)*(n(:)+1).^(exponent3y(2));
                        % find second zero
                        x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
                    end
                otherwise
                    error('Code should never get here.');
            end
            if k >= 4
                for iRoot = 4:k
                    % guesses for remaining roots x[k] = rootSpacing + x[k-1]
                    x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = 2*x(((iRoot-2)*nOrdersPerRoot+1:(iRoot-1)*nOrdersPerRoot)')- x(((iRoot-3)*nOrdersPerRoot+1:(iRoot-2)*nOrdersPerRoot)');
                    % find the remaining zeros
                    x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = arrayfun(@(n, x0) obj.findzero(n, k, x0, kind), n, x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)'));
                end
            end
            end
            function x=findzero(obj,n,k,x0,kind)
            % Uses Halley's method to find a zero given the starting point x0
            % http://en.wikipedia.org/wiki/Halley's_method
            ITERATIONS_MAX = 100;       % Maximum number of iteration
            TOLERANCE_RELATIVE = 1e4;   % 16-4 = 12 significant digits
            % Setup loop
            error = 1;
            loopCount = 0;
            x = 1; % Initialization value only.  It is is not used.
            % Begin loop
            while abs(error)>eps(x)*TOLERANCE_RELATIVE && loopCount<ITERATIONS_MAX

                switch kind
                    case 1
                        a = besselj(n,x0);
                        b = besselj((n+1),x0);
                    case 2
                        a = bessely(n,x0);
                        b = bessely((n+1),x0);
                end

                xSquared = x0*x0;

                error = 2*a*x0*(n*a-b*x0)/(2*b*b*xSquared-a*b*x0*(4*n+1)+(n*(n+1)+xSquared)*a*a);

                % Prepare for next loop
                x=x0-error;
                x0=x;
                loopCount=loopCount+1;

            end
            % Handle maximum iterations
            if loopCount>ITERATIONS_MAX-1
                warning('Failed to converge to within relative tolerance of %e for n=%f and k=%d in %d iterations', eps(x)*TOLERANCE_RELATIVE, n, k, ITERATIONS_MAX);
            end
            end

        end
        
 end 
