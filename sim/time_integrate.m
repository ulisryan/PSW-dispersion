function velocity = time_integrate(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend, use_gpu)
% TIME_INTEGRATE  Velocity Verlet 시간 적분.
%
% use_gpu 인자(선택):
%   true  — GPU 사용 (단독 실행 시)
%   false — CPU 사용 (parfor worker 내부에서 호출 시 권장)
%   생략  — GPU 자동 감지

    if nargin < 17
        use_gpu = gpuDeviceCount('available') > 0;
    end

    if use_gpu
        %% GPU 경로: VRAM 크기에 맞게 자동 청킹
        g_dev = gpuDevice;
        nI = size(aligned_cost, 2);
        bytes_per_particle = nI * 8 + 250;
        gpu_batch_size = floor(g_dev.AvailableMemory * 0.8 / bytes_per_particle);
        gpu_batch_size = min(gpu_batch_size, MAX);

        n_batches = ceil(MAX / gpu_batch_size);
        fprintf('  GPU: %s | %.2f GB 가용 | %d개씩 %d배치\n', ...
            g_dev.Name, g_dev.AvailableMemory/1e9, gpu_batch_size, n_batches);

        tic;
        for b = 1:n_batches
            idx   = ((b-1)*gpu_batch_size + 1) : min(b*gpu_batch_size, MAX);
            b_sz  = numel(idx);
            velocity(idx,:) = gpu_batch( ...
                position(idx,:), velocity(idx,:), Ij(idx), Ips(idx,:), ...
                I1_vec(idx), I2_vec(idx), aligned_cost(idx,:), last_data(idx,:), ...
                lambda, w0z, w0y, tau, b_sz, tstep, t_ini, tend);
        end
        fprintf('  완료 %.1f초\n', toc);

    else
        %% CPU 경로: BLAS 멀티스레드 (parfor worker 내 또는 GPU 없는 환경)
        t = t_ini;
        while t < tend
            position = position + velocity .* tstep;
            [~, last_data, a, ~] = axay_opt_partial( ...
                position, t + Ij, I1_vec, I2_vec, ...
                aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX);
            velocity = velocity + a .* tstep;
            t = t + tstep;
        end
    end
end


function velocity = gpu_batch(position, velocity, Ij, Ips, I1_vec, I2_vec, ...
        aligned_cost, last_data, lambda, w0z, w0y, tau, MAX, tstep, t_ini, tend)
    position     = gpuArray(position);
    velocity     = gpuArray(velocity);
    aligned_cost = gpuArray(aligned_cost);
    I1_vec       = gpuArray(I1_vec);
    I2_vec       = gpuArray(I2_vec);
    Ij           = gpuArray(Ij);
    Ips          = gpuArray(Ips);
    last_data    = gpuArray(last_data);

    t = t_ini;
    while t < tend
        position = position + velocity .* tstep;
        [~, last_data, a, ~] = axay_opt_partial( ...
            position, t + Ij, I1_vec, I2_vec, ...
            aligned_cost, last_data, Ips, lambda, w0z, w0y, tau, MAX);
        velocity = velocity + a .* tstep;
        t = t + tstep;
    end
    velocity = gather(velocity);
end
