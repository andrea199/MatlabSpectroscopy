classdef ADEModule < handle
    %Class with routines to run ADE with Taylor Dispersion
    %c_t + u*c_x = ((R/L)^2* Pe/48 *u^2 + 1/Pe)*c_xx
    %var - Up, fun - Low
    
    properties
        RoLScale = 0.1; % Appropriate radius-to-length scale in the problem -> updated from GeometryParserModule
        PecletNum = 1;  % Peclet numer -> updated from GeometryParserModule
        UScale = 100; % Velocity scale -> updatseed from GeoemtryParserModule
        Ngrid = 200;  % Number of grid points -> updated from GUI
        DomainLength = 10; %Domain length -> updated from GeometryParserModule
        RealDomainIndex = 1; % Index marking the end real domain (to accurately model the outflow BC)
        DeltaT = 0.0001; %Time step -> updated from GeometryParserModule (From pump make n model)
        DcoeffFn; % Function handle to hold the defintion of the diffusion coefficient
        Dcoeff = []; % Peclet number with Taylor Diffusion
        UField = []; %Depth averaged velocity field -> updated from GeometryParserModule
        Solver; % Structure for holding the solver properties -> partly updated from GeometryParserModule
    end
    
   methods % Constructor
        function obj = ADEModule(SolverSelect)
          % Handle optional arguments
          if nargin<1
              SolverSelect = 'CN';
          end 
          obj.Solver.SolverSelect = SolverSelect;      
       end     
   end 
     
   methods % Core functions:
        function stepSpectralSolver(obj)
            % Evolves the solution by one time step
            obj.Solver.q=obj.Solver.M*obj.Solver.q - obj.Solver.G\(obj.Solver.Alpha*obj.DeltaT); 
        end
        function stepCrankNicolson(obj)
            %Evaluates the solution by one time step
            obj.Solver.q =(obj.Solver.M\obj.Solver.N)*obj.Solver.q; 
        end
        function stepBackwardEuler(obj)
            %Evaluates the solution by one time step
            obj.Solver.q =(obj.Solver.BE.M\obj.Solver.BE.N)*obj.Solver.q; 
        end     
   end
   
   methods %Supporting functions: Spectral Solver
       function initializeCrankNicolson(obj)
           % Initialize the CrankNicolson solver  wit h the Neumann boundary conditions at
           % the delivery site, and part of the Dirchlet boundary condition
           % at the source
           
            
           obj.UField = obj.UScale;
           % Assemble Matrices
           DeltaZ= obj.DomainLength/obj.Ngrid;
           Alpha = obj.UField.*obj.DeltaT/(4*DeltaZ);
           Beta  = obj.Dcoeff.*obj.DeltaT/(2*DeltaZ^2); 
           
           %Sparse matrices (old version)
%           M = spdiags([-Alpha-Beta, 1+2*Beta, Alpha-Beta],-1:1,obj.Ngrid,obj.Ngrid);
%           N = spdiags([Alpha+Beta, 1-2*Beta, -Alpha+Beta],-1:1,obj.Ngrid,obj.Ngrid);
%             
           %Sparse matrices (Corrected to match python performance)
            M = spdiags([Alpha-Beta, 1+2*Beta, -Alpha-Beta],-1:1,obj.Ngrid,obj.Ngrid);
            N = spdiags([-Alpha+Beta, 1-2*Beta, Alpha+Beta],-1:1,obj.Ngrid,obj.Ngrid);
            M = M';
            N = N';
           
           %Bounday conditions
           M(1,1) = 1;
           M(1,2) = 0;
           M(obj.Ngrid,obj.Ngrid) = 1;
           M(obj.Ngrid,obj.Ngrid-1) = -1;
           
           N(1,1) = 1;
           N(1,2) = 0;
           N(obj.Ngrid,obj.Ngrid) = 0;
           N(obj.Ngrid,obj.Ngrid-1) = 0;
           
           % Assign Matrices
           obj.Solver.M = M;
           
           % Updated definition of N so that U and D are functions of time
           % Change - 11/3/2022
           if ~isfield(obj.Solver, 'Nprev')
           obj.Solver.N = N;
           else
           obj.Solver.N = obj.Solver.Nprev;               
           end
           obj.Solver.Nprev=N;
           
           obj.Solver.x = linspace(0,obj.DomainLength,obj.Ngrid);
           % Initial and Bounday conditions
            obj.Solver.qint = zeros(obj.Ngrid,1); 
            obj.Solver.qint(1)=1; 
            obj.Solver.q = obj.Solver.qint;    
            
       end      
       function initializeBackwardEuler(obj)
           % Initialize the BackwardEuler solver  with the Neumann boundary conditions at
           % the delivery site, and part of the Dirchlet boundary condition
           % at the source. 
           % Purpose: CN will be supplemented with BE to control
           % oscillations.
           
           obj.UField = obj.UScale;
           % Assemble Matrices
           DeltaZ= obj.DomainLength/obj.Ngrid;
           Alpha = obj.UField.*obj.DeltaT/(2*DeltaZ);
           Beta  = obj.Dcoeff.*obj.DeltaT/(DeltaZ^2); 
           
           %Sparse matrixs
           M = spdiags([-Alpha-Beta, 1+2*Beta, Alpha-Beta],-1:1,obj.Ngrid,obj.Ngrid);     
           N = speye(size(M));
           %Bounday conditions
           M(1,1) = 1;
           M(1,2) = 0;
           M(obj.Ngrid,obj.Ngrid) = 1;
           M(obj.Ngrid,obj.Ngrid-1) = -1;  
           N(obj.Ngrid,obj.Ngrid) = 0;
           %Assign Matrices
           obj.Solver.BE.M = M;   
           obj.Solver.BE.N = N;
       end
   end
   
   
   methods % Derived methods for higher order calculations
       % Intended functions: 
       % (i) Cumulative multi drug solver (injection port closest to
       % patient will have the contribution from all the drugs, the one before
       % will have one less drug and so on.). For implementation each port
       % will have it's own solver with number of qd = no. of drugs
       % immediately down stream)
       % (ii) Scaling the variable such that at x=0, c = 1. This will
       % simplify the definition of the objective function used for the optimization.
       % Decided to define a new class.
     
       
   end
   
   methods %Supporting functions: Spectral Solver
       function initializeSpectralSolver(obj)
           % Initialize the spectral solver  with the Neumann boundary conditions at
           % the delivery site, and part of the Dirchlet boundary condition
           % at the source
           
           % Generate differentiation matrices
           obj.Solver.scale = -2/obj.DomainLength;
           [obj.Solver.x, obj.Solver.DM] = obj.chebdif(obj.Ngrid,2); % We are solving a second order ODE
           obj.Solver.dx=obj.Solver.DM(:,:,1)*obj.Solver.scale;
           obj.Solver.dxx=obj.Solver.DM(:,:,2)*obj.Solver.scale^2;
           obj.Solver.x=(obj.Solver.x-1)/obj.Solver.scale; 
           obj.Solver.I=eye(obj.Ngrid);
           
           % System Matrices
            obj.Solver.E=obj.Solver.I;
            obj.Solver.A=-obj.UField.*obj.Solver.dx+obj.Dcoeff.*obj.Solver.dxx; %The depth averaged velocity field can vary in space.

            % Boundary conditions
            obj.Solver.E([1 obj.Ngrid],:)=0;
            obj.Solver.A(1,:)=obj.Solver.I(1,:);
            obj.Solver.A(obj.Ngrid,:)=obj.Solver.dx(obj.Ngrid,:);
            obj.Solver.Alpha = obj.Solver.I(:,1); %The homogeneous boundary condition
            % March in time matrix 
            obj.Solver.M=(obj.Solver.E-obj.Solver.A*obj.DeltaT/2)\(obj.Solver.E+obj.Solver.A*obj.DeltaT/2) ;
            obj.Solver.G = (obj.Solver.E-obj.Solver.A*obj.DeltaT/2);
            % Initial conditions
            obj.Solver.qint = zeros(obj.Ngrid,1);  
            obj.Solver.q = obj.Solver.qint;
       end    
       function [x, DM] = chebdif(~, N, M)

        %  The function DM =  chebdif(N,M) computes the differentiation 
        %  matrices D1, D2, ..., DM on Chebyshev nodes. 
        % 
        %  Input:
        %  N:        Size of differentiation matrix.        
        %  M:        Number of derivatives required (integer).
        %  Note:     0 < M <= N-1.
        %
        %  Output:
        %  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
        %
        %  The code implements two strategies for enhanced 
        %  accuracy suggested by W. Don and S. Solomonoff in 
        %  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
        %  The two strategies are (a) the use of trigonometric 
        %  identities to avoid the computation of differences 
        %  x(k)-x(j) and (b) the use of the "flipping trick"
        %  which is necessary since sin t can be computed to high
        %  relative precision when t is small whereas sin (pi-t) cannot.
    


     I = eye(N);                          % Identity matrix.     
     L = logical(I);                      % Logical identity matrix.

    n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

     k = [0:N-1]';                        % Compute theta vector.
    th = k*pi/(N-1);

     x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

     T = repmat(th/2,1,N);                
     DX = 2*sin(T'+T).*sin(T'-T);          % Trigonometric identity. 
     DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];   % Flipping trick. 
     DX(L) = ones(N,1);                       % Put 1's on the main diagonal of DX.

     C = toeplitz((-1).^k);               % C is the matrix with 
     C(1,:) = C(1,:)*2; C(N,:) = C(N,:)*2;     % entries c(k)/c(j)
     C(:,1) = C(:,1)/2; C(:,N) = C(:,N)/2;

     Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))  
     Z(L) = zeros(N,1);                      % with zeros on the diagonal.

     D = eye(N);                          % D contains diff. matrices.
                                          
    for ell = 1:M
              D = ell*Z.*(C.*repmat(diag(D),1,N) - D); % Off-diagonals
           D(L) = -sum(D');                            % Correct main diagonal of D
    DM(:,:,ell) = D;                                   % Store current D in DM
    end
       end
   end
end