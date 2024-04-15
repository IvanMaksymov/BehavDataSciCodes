###############################################################################
# This is an Octave script that solves the Schrodinger equation
# in two-dimensional space using the Crank-Nicolson method.
# The script draws on the Python code available at
# https://artmenlope.github.io/solving-the-2d-schrodinger-equation-using-the-crank-nicolson-method/
# Matlab users will need to substitutre the comment characters # for %.
# All parameters are chosen so that the computation takes a reasonable for
# demonstration purposes time. The users are invited to change the domain size
# L and adjust the number of # time steps Nt to obtain high-resolution images.
###############################################################################

clc;close all;clear all;

fnt_size = 14; #font size for the figures

# Setting the dimensions of the computational domain
L = 8;   dy = 0.05;   dt = dy^2/4;
N = round(L/dy) + 1;   Nt = 500;
rx = -dt/(2*1i*dy^2);   ry = -dt/(2*1i*dy^2);

x0 = L/5;   y0 = L/2;  #x0 = L/2;   y0 = L/2;
x01 = 4*L/5;   y01 = L/2;
x = linspace(0,L,N);   y = linspace(0,L,N);   [X, Y] = meshgrid(x,y);

# Setting the incodent Gaussian pulse (uncomment the figure to view its shape)
psi0 = zeros(N,N);   sigma = 0.5;   sigma1 = 0.5;   k=15*pi;   k1=15*pi;
#change parameters sigma and k to modify the original Gaussian shapes of the pulses

#two Gaussian pulses
psi0 = exp(-1/2*((X-x0).^2 + (Y-y0).^2)/sigma^2).*exp(1i*k*(X-x0)) + ...
       exp(-1/2*((X-x01).^2 + (Y-y01).^2)/sigma1^2).*exp(-1i*k1*(X-x01));

#one Gaussian pulse
#psi0 = exp(-1/2*((X-x0).^2 + (Y-y0).^2)/sigma^2).*exp(1i*k*(X-x0));

#figure(10);subplot(1,1,1);imagesc(abs(psi0));axis square;

# Double slit structure
w = 0.6;   s = 0.8;   a = 0.2;

j0 = round(1/(2*dy)*(L-w));   j1 = round(1/(2*dy)*(L+w));

i0 = round(1/(2*dy)*(L+s) + a/dy);   i1 = round(1/(2*dy)*(L+s));
i2 = round(1/(2*dy)*(L-s));   i3 = round(1/(2*dy)*(L-s) - a/dy);

v0 = 0;#0 -- no double-slit screen, 200 -- typical value, 2000 and higher -- "impermeable" screen
v = zeros(N,N);   v(1:i3,j0:j1) = v0;   v(i2:i1,j0:j1) = v0;   v(i0:end,j0:j1) = v0;


# Uncomment to inspect the double-slit structure
# figure(1);subplot(1,2,1);imagesc(v);axis square;
# rectangle ("Position", [73.5, 0.5, 13, 68.5]);
# rectangle ("Position", [73.5, 71.5, 13, 17]);
# rectangle ("Position", [73.5, 91.5, 13, 69.5]);
# drawnow;


Ni = N*N;   A = zeros(Ni,Ni);   M = zeros(Ni,Ni);

# Populating the elements of the main matrix
k = 0;
for j = 1:N
  for i = 1:N

      k = k + 1;

      A(k,k) = 1 + 2*rx + 2*ry + 1i*dt/2*v(i,j);
      M(k,k) = 1 - 2*rx - 2*ry - 1i*dt/2*v(i,j);

      if k != 1
        A(k-1,k) = -rx;   M(k-1,k) = rx;
      end;

      if k != Ni
        A(k+1,k) = -rx;   M(k+1,k) = rx;
      end;

      if k > N
        A(k-N,k) = -ry;   M(k-N,k) = ry;
        A(k,k-N) = -ry;   M(k,k-N) = ry;
      end;

  end;
end;

# Uncomment to inspect the main matrix
#figure(1);subplot(1,2,2);spy(abs(A),'ro',1);axis square;drawnow;

psi = zeros(N,N,Nt);   psi = psi0;
# Boundary conditions
psi(1,1:end) = 0.0;   psi(1:end,1) = 0.0;   psi(end,1:end) = 0.0;   psi(1:end,end) = 0.0;

psi_vect = reshape(psi,Ni,1);   Asp = sparse(A);   Msp = sparse(M);

# Solution in the time domain
for i = 1:Nt

    b = mtimes(Msp,psi_vect);   psi_vect = Asp\b;#linsolve(A,b);
    psi(1,1:end) = 0.0;   psi(1:end,1) = 0.0;   psi(end,1:end) = 0.0;   psi(1:end,end) = 0.0;
    psi(:,:,i) = reshape(psi_vect,N,N);

end;

# Record videos
# Uncomment this part of the code to enable video recording
# The user will need to download and install open FFMPEG software (see https://ffmpeg.org/)
# The following section of the code works in Linux/Mac.
# Some changes may needed to be by Windows users

out_dir = "temp_img";   mkdir (out_dir);

for i = 1:Nt

  figure(1000);subplot(1,1,1);hold on;imagesc(abs(psi(:,:,i)).^2);axis square;
  caxis([0 1.0]);   xlim([-10 170]);   ylim([-10 170]);
  #rectangle ("Position", [73.5, 0.5, 13, 68.5]);
  #rectangle ("Position", [73.5, 71.5, 13, 17]);
  #rectangle ("Position", [73.5, 91.5, 13, 69.5]);
  set(gca,'fontsize',fnt_size);   set(gca,'linewidth',2);
  xlabel('x-coordinate (arb. units)');   ylabel('y-coordinate (arb. units)');
  colorbar;   box on;

  fname = fullfile (out_dir, sprintf ("img%03i.png", i));
  imwrite (getframe (gcf).cdata, fname);

  drawnow;

end;

cmd = sprintf ("ffmpeg -framerate 20 -i ./%s/img%%03d.png -vf scale=1080:-1 Example1.mp4", out_dir);
system (cmd)

STOP


