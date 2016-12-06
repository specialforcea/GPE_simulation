function [Dipolar_interaction_nonlinearity] = Dipolar_interaction3d(Phi, FFTX, FFTY, FFTZ, Dipolar_direction, d)
% Computing the norm of the cross product between the frequency vector and the dipolar direction
Cross_norm = sqrt((FFTY*Dipolar_direction(3)-FFTZ*Dipolar_direction(2)).^2+(FFTZ*Dipolar_direction(1)-FFTX*Dipolar_direction(3)).^2+(FFTX*Dipolar_direction(2)-FFTY*Dipolar_direction(1)).^2);
% Computing the norm of the scalar product between the frequency vector and the dipolar direction
Scalar_prod = FFTX*Dipolar_direction(1)+FFTY*Dipolar_direction(2)+FFTZ*Dipolar_direction(3);
% Computing the angle between the frequency vector and the dipolar
% direction
Angle = atan2(Cross_norm,Scalar_prod);
% Computing the FFT of the square of the module of phi
NLFFT = fftn(abs(Phi).^2);
% Computing the FFT potential
V = d^2*(4/3)*pi*(3*cos(Angle).^2-1);
% Computing the dipolar interaction
Dipolar_interaction_nonlinearity = ifftn(V.*NLFFT);