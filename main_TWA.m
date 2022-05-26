%MAIN FILE FOR THE 2D HEISENBERG MODEL USING TWA:

%Description: code to compute the dynamics of spins in a 2D cubic lattice
%using the Truncated Wigner Approximation and a spin spiral as initial
%condition. We use the RK4 method as integrator. 
%Output: we compute as an example the equal-time Fourier transform of the
%spin components averaged radially in k-space.


%Model parameters (see manuscript for details):
L = 100;%Linear size
S = 10;%Spin number
qx = 8*(2*pi()/L);%Wavevector of the spin spiral
qy = 0;%Wavevector of the spin spiral
theta = pi()/2;%Angle of the spin spiral 
Delta = 0.0;%Spin anisotropy [Jz = J(1+Delta)]

%Simulation parameters:
N_samples = 1;
psi_0 = sin(theta);

%Time range:
tau = 1/(psi_0^2*(2-cos(qx)-cos(qy)));%define characteristic timescale
T = tau;%Time window between compute times
RPT = 10;%Number of 
N_steps = 1000*(T/tau);%Number of steps for the RK integration
nsigma = 2;

Ck_xx = zeros(length(kr_range),RPT+1);
Ck_zz = zeros(length(kr_range),RPT+1);

fprintf('---Start Sampling---\n\n')
for sample = 1 : N_samples

    %Create sample:
    psi_x = zeros(L,L);
    psi_y = zeros(L,L);
    psi_z = zeros(L,L);
    [psi_x,psi_y,psi_z] = create_spiralgaussian(L,qx,qy,theta,S);
    
    %Compute:
    [Cii,kr_range] = compute_FT(psi_x,nsigma);
    Ck_xx(:,1) = Ck_xx(:,1) + Cii;
    [Cii,kr_range] = compute_FT(psi_z,nsigma);
    Ck_zz(:,1) = Ck_zz(:,1) + Cii;
    
    %run
    for repeat = 1 : RPT
        tic
        [psi_x,psi_y,psi_z] = run_twark4(psi_x,psi_y,psi_z,Delta,T,N_steps);
        
        %Compute:       
        [Cii,kr_range] = compute_FT(psi_x,nsigma);
        Ck_xx(:,repeat+1) = Ck_xx(:,repeat+1) + Cii;
        [Cii,kr_range] = compute_FT(psi_z,nsigma);
        Ck_zz(:,repeat+1) = Ck_zz(:,repeat+1) + Cii;

        toc
        fprintf('Sampling percentage: %6.2f.\n',sample/N_samples*100)
        fprintf('Step %6.2f.\n\n',repeat/RPT*100)
    end
end

%filename = %Filename
%save(filename,'Ck_xx','Ck_zz')