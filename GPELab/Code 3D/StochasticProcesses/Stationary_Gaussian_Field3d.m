function [B] = Stationary_Gaussian_Field3d(Geometry3D,varargin)
%% Analysis of the inputs
Analyse_Var = inputParser; % Creating the parser
Analyse_Var.addOptional('Cov_fun', @(X,Y,Z) exp(-X.^2-Y.^2-Z.^2)); % Required input 'Cov_fun' corresponding to a covariance function

%% Computing the geometry for the FFT
FFTGeometry3D = FFTGeometry3D_Var3d(Geometry3D);

%% Parsing inputs and storing inputs
% Parsing inputs
Analyse_Var.parse(varargin{:}); % Analysing the inputs
Cov_fun = Analyse_Var.Results.Cov_fun; %Storing the 'Method' input

% Computing a random normal matrix
RandN = randn(FFTGeometry3D.Ny,FFTGeometry3D.Nx,FFTGeometry3D.Nz);

% Computing the gaussian field using the spectral representation
Cov = Cov_fun(FFTGeometry3D.X,FFTGeometry3D.Y,FFTGeometry3D.Z);
FFTRandN = fftn(RandN);
FFTCov = fftn(fftshift(Cov));
B_tmp = real(ifftn(sqrt(FFTCov).*FFTRandN));

% Adding periodic boundaries
B = zeros(FFTGeometry3D.Ny+1,FFTGeometry3D.Nx+1,FFTGeometry3D.Nz+1);
B(1:FFTGeometry3D.Ny,1:FFTGeometry3D.Nx,1:FFTGeometry3D.Nz) = B_tmp;
B(FFTGeometry3D.Ny+1,1:FFTGeometry3D.Nx,1:FFTGeometry3D.Nz) = B_tmp(1,:,:);
B(1:FFTGeometry3D.Ny,FFTGeometry3D.Nx+1,1:FFTGeometry3D.Nz) = B_tmp(:,1,:);
B(1:FFTGeometry3D.Ny,1:FFTGeometry3D.Nx,FFTGeometry3D.Nz+1) = B_tmp(:,:,1);
B(FFTGeometry3D.Ny+1,FFTGeometry3D.Nx+1,FFTGeometry3D.Nz+1) = B_tmp(1,1,1);
