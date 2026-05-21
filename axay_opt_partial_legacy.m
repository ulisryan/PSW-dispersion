function [cos2t,last_data,acceleartion,temp] = axay_opt_partial(position, t, I1,I2, aligned_value, last_data, IR_position_stability,lambda,w0z,w0y,tau,MAX)
    % Constants
    upsilon = 8.854187e-12;
    speed_light = 2.998e8;
    % lambda = 1064e-9;    % Laser wavelength
    % w0 = 25e-6;        % Waist radius of laser
    % tau = 7.54e-9;       % Laser pulse width
    q = 2 .* pi .* position(:,1) ./ lambda;
    n = 2 * pi ./ lambda;
    %mm = 0;
    acceleartion=zeros(MAX,3);
    a_list=zeros(MAX,3);
    % Calculate laser intensity I
    I = (I1 + I2 + 2 .* sqrt(I1 .* I2) .* cos(2.*q)) .* exp((-2 .* (position(:,3) + IR_position_stability(:,1)).^2) ./ w0z^2) .* exp(-2 .* (position(:,2)+IR_position_stability(:,2)).^2 ./ w0y^2); % waist radius 이기 때문에 2를 곱해주는것이 맞다.
    % Interpolate laser intensity for aligned value
    interpolating_laser_intensity = I .* exp(-4 .* log(2.0) .* (t.^2) ./ (tau.^2)) .* (1e-13);
    interpolate_I = floor(interpolating_laser_intensity); % Interpolate with floor value (ex. 2.6333e14 -> 2.6e14) 

     row_indices=(1:MAX)';
     column_indices=interpolate_I+1;
     key=sub2ind(size(aligned_value), row_indices, column_indices);
     
    cos2t = aligned_value(key) + (interpolating_laser_intensity - interpolate_I) .* (aligned_value(key + MAX) - aligned_value(key));
    
    % Calculate alpha
    alpa = ((16.8 - 6.2) * 1e-40 * cos2t + 6.2 * 1e-40);
    mm = 0.5 * I / upsilon / speed_light;

    %if last_data(:,2) ~= 0
    % last_data
        logic=last_data(:,2)~=0;
        % logic
        da_dI = (last_data(:,1) - alpa) ./ (last_data(:,2) - I .* exp(-4 * log(2.0) .* t.^2 / tau^2));
    % da_dI
        da_dI = da_dI.*logic;
    %end
   % da_dI
    % Calculate acceleration components (I~=0)
a_list(:,1) = -mm .* 2.* n .* ( 2 .* sqrt(I1.*I2) .* (sin(2.*q)) ./ (I1+I2+2.*sqrt(I1.*I2).*cos(2.*q)) ) ./ (76.14e-3 ./ 6.022e23) .* exp(-4 .* log(2.0) .* t .* t / tau / tau) .* (alpa + I .* da_dI .* exp(-4 .* log(2.0) .* t .* t / tau / tau));
a_list(:,2) = -4 .* mm .* position(:,2) / (76.14e-3 / 6.022e23) .* exp(-4 .* log(2.0) .* t .* t / tau / tau) / w0y / w0y .* (alpa + I .* da_dI .* exp(-4 .* log(2.0) .* t .* t / tau / tau));
a_list(:,3) = -4 .* mm .* position(:,3) / (76.14e-3 / 6.022e23) .* exp(-4 .* log(2.0) .* t .* t / tau / tau) / w0z / w0z .* (alpa + I .* da_dI .* exp(-4 .* log(2.0) .* t .* t / tau / tau));
% a_list
acceleartion(I~=0,1) = a_list(I~=0,1);
acceleartion(I~=0,2) = a_list(I~=0,2);
acceleartion(I~=0,3) = a_list(I~=0,3);

acceleartion(I==0,1) =0;
acceleartion(I==0,2) =0;
acceleartion(I==0,3) =0;
% acceleartion
    % Update last_data
    last_data(:,1) = alpa;
    last_data(:,2) = I .* exp(-4 .* log(2.0) .* t.^2 / tau.^2);

    %임시로 함수 확인을 위해 받는 데이터
    temp=[I,interpolating_laser_intensity,interpolate_I,row_indices,column_indices,key,cos2t];
end