fid = fopen(jm_data, 'r');

alpha_data=fscanf(fid, '%lf', [4, 1224]);
fclose(fid);
alpha_data=alpha_data';
input_j=alpha_data(:,1);
input_m=alpha_data(:,2);
alpha=alpha_data(:,3);
proba=alpha_data(:,4);
proba(2:end)=proba(2:end)-proba(1:end-1);

%random distribution for j m
random_index=randsample(1:1224,MAX,true,proba);
particle(:,10)=input_j(random_index);
particle(:,11)=abs(input_m(random_index));

%load polarizability data
aligned_cost=zeros(MAX,2001);
load('linearly.mat')

%postion and velocity
particle(:,1)=x_position_fwhm/2.355 *randn(MAX,1); 
% particle(:,1)=-500+randi([0 999],MAX,1); % -500 ~ 500

particle(:,2)=Dye_sigma*randn(MAX,1); 
particle(:,4)=vx_caused_tilted+x_init_vel_sigma*randn(MAX,1);
particle(:,5)=y_init_vel_sigma*randn(MAX,1);
particle(:,6)=z_init_vel_off+z_init_vel_sigma*randn(MAX,1);

laser(:,1)=Dye_sigma*randn(MAX,1);
laser(:,2)=(Dye_tau)*randn(MAX,1)/sqrt(3);
laser(:,3)= (-Ij_distribution+randi([0 Ij_distribution*2],MAX,1))* 1e-12;
laser(:,4)= Ips_z.*randn(MAX,1)* 1e-9;
laser(:,5)= Ips_y.*randn(MAX,1)* 1e-9;

particle(:,3) = (30e-9 + laser(:,2)) .* z_init_vel_off - (60e-9 + laser(:,2)) .* particle(:,6) + laser(:,1);

last_data=zeros(MAX,2);
position=particle(:,1:3);
velocity=particle(:,4:6);
acceleration=zeros(MAX,3);
Ij=laser(:,3);
Ips=[laser(:,4),laser(:,5)];

% load linearly
for index=1:MAX
aligned_cost(index,:)=linearly(particle(index,10)+1,particle(index,11)+1,:);
end

max_I = 2.0 / (pi * w0z * w0y) * sqrt(4.0 * log(2.0) / (pi * tau * tau));
% max_I = max_I * 0.959 *0.948*1.0e-3;
max_I = max_I * atten1 *atten2*1.0e-3;
I1=max_I .* I1_list(I_ind);
I2=max_I .* I2_list(I_ind);

tcount=0;
t=t_ini;

tic
while t<tend

    position=position+velocity.*tstep;
    [~,last_data, a,temp] = axay_opt_partial(position, t + Ij, I1,I2, aligned_cost, last_data, Ips,lambda,w0z,w0y,tau,MAX);
    velocity=velocity+a.*tstep;
    t = t + tstep;
    tcount=tcount+1;

end
toc
