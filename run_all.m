% run_all.m — PSW Dispersion Simulation
%
% 병렬 전략: parfor를 sub-batch 단위로 적용
%   MAX=200,000 × 8조건 → sub-batch 20,000 × 80태스크
%   8 workers (메모리 대역폭 포화점, 실측 최적)
%   실측 시간: ~181초 (기존 순차 ~960초 대비 5.3배 빠름)

clear; clc;
tic

root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir, fullfile(root_dir, 'sim'));

%% ── 0. Parallel pool 설정 (MaxNumWorkers 강제 확장) ──────────────────────────
n_workers = 8;   % 실측 최적값 (8=181s, 12=183s, 16=187s, 24=208s)

% 기존 풀이 worker 수가 부족하면 재시작
p = gcp('nocreate');
if ~isempty(p) && p.NumWorkers ~= n_workers
    delete(p);
    p = [];
end

if isempty(p)
    cluster = parcluster('local');
    if cluster.NumWorkers ~= n_workers
        cluster.NumWorkers = n_workers;  % 프로파일 MaxNumWorkers 고정
        saveProfile(cluster);            % 변경사항 저장
    end
    parpool(cluster, n_workers);
    fprintf('Parallel pool 시작: %d workers\n', n_workers);
end

%% ── 1. 상수 로드 ─────────────────────────────────────────────────────────────
load(fullfile(root_dir, 'STN_constants.mat'));
MAX = 2000000;   % 조건당 총 입자 수 (STN_constants의 값 overwrite)

%% ── 2. 이미지 처리 파라미터 ──────────────────────────────────────────────────
rep_std       = 4.0;
rep_vel       = 2.7;
vel_pix       = 2.28;
atten_factorx = 1 / 0.9278;
atten_factory = 1 / 0.9453;

%% ── 3. 양자상태 데이터 파일 경로 ─────────────────────────────────────────────
jm_data_path = fullfile(root_dir, 'data', jm_data);
if ~isfile(jm_data_path)
    error('파일 없음: %s', jm_data_path);
end

%% ── 4. linearly.mat 사전 로드 (parfor broadcast 변수) ────────────────────────
fprintf('[1/2] linearly.mat 로드 중...\n');
load(fullfile(root_dir, 'linearly.mat'));
[nj_lin, nm_lin, ~] = size(linearly);
linearly_2d = reshape(linearly, nj_lin*nm_lin, []);  % (nj*nm)×nI — 모든 worker에 공유
clear linearly;

%% ── 5. 레이저 세기 사전 계산 ─────────────────────────────────────────────────
max_I   = 2.0 / (pi * w0z * w0y) * sqrt(4.0 * log(2) / (pi * tau^2));
max_I   = max_I * atten1 * atten2 * 1e-3;
I1_val  = max_I * I1_list;       % 스칼라
I2_vals = max_I * I2_list(:);   % n_groups×1

n_groups   = numel(I2_list);
batch_size = 20000;                          % sub-batch당 입자 수 (GPU VRAM 또는 캐시 기준)
n_sub      = ceil(MAX / batch_size);         % 조건당 sub-batch 수 (예: 10)
n_tasks    = n_groups * n_sub;               % 전체 parfor 태스크 수 (예: 80)

fprintf('[2/2] Time integration 시작...\n');
fprintf('      조건: %d × sub-batch: %d × 입자: %d = 총 태스크 %d\n', ...
    n_groups, n_sub, batch_size, n_tasks);
fprintf('      24 logical core → ceil(%d/24)=%d 라운드 예상\n', ...
    n_tasks, ceil(n_tasks/24));

%% ── 6. parfor: sub-batch 단위 병렬 시뮬레이션 ───────────────────────────────
% task 번호 → (조건 g, sub-batch s) 매핑:
%   task = (g-1)*n_sub + s
%   → g = ceil(task/n_sub),  s = mod(task-1, n_sub)+1

results = cell(n_tasks, 1);

parfor task = 1:n_tasks

    g_idx = ceil(task / n_sub);   % 어떤 I2 조건인지 알아내기. 결국 총 task에서 한 세기당의 몇 sub batch인지를 나누어 주면 몇번째 세기인지 알 수 있음.

    %% 입자 초기화 (sub-batch 크기만큼)
    [particle_t, laser_t] = init_particles(jm_data_path, batch_size, ...
        x_position_fwhm, Dye_sigma, Dye_tau, ...         %#ok<PFBNS>
        vx_caused_tilted, x_init_vel_sigma, y_init_vel_sigma, ...
        z_init_vel_off,   z_init_vel_sigma, ...
        Ij_distribution,  Ips_z, Ips_y);

    %% aligned_cost 벡터화 생성
    j_idx_t = particle_t(:,10) + 1;
    m_idx_t = particle_t(:,11) + 1;
    aligned_cost_t = linearly_2d(sub2ind([nj_lin, nm_lin], j_idx_t, m_idx_t), :);

    %% 레이저 세기 벡터
    I1_t = repmat(I1_val,         batch_size, 1);
    I2_t = repmat(I2_vals(g_idx), batch_size, 1);

    %% CPU 시간 적분
    results{task} = time_integrate( ...
        particle_t(:,1:3), particle_t(:,4:6), ...
        laser_t(:,3),      laser_t(:,4:5), ...
        I1_t, I2_t, aligned_cost_t, zeros(batch_size, 2), ...
        lambda, w0z, w0y, tau, batch_size, tstep, t_ini, tend, false);

end

%% ── 7. sub-batch 결과를 조건별로 재조립 ────────────────────────────────────
Cell_velocity = cell(n_groups, 1);
for g = 1:n_groups
    task_idx = (g-1)*n_sub + (1:n_sub);           % 이 조건에 속한 task 번호들
    Cell_velocity{g} = vertcat(results{task_idx}); % 200,000×3로 합치기
end

%% ── 8. 이미지 처리 (parfor) ──────────────────────────────────────────────────
% fprintf('[3/3] Image processing (parfor)...\n');
% [img_matrix, sumx, sumy] = image_processing(Cell_velocity, I2_list, ...
%     rep_std, rep_vel, vel_pix, atten_factorx, atten_factory);

%% ── 9. 결과 저장 ─────────────────────────────────────────────────────────────
% 시뮬레이션 제어
constants.n_workers     = n_workers;
constants.MAX           = MAX;
constants.batch_size    = batch_size;
% 이미지 처리
constants.rep_std       = rep_std;
constants.rep_vel       = rep_vel;
constants.vel_pix       = vel_pix;
constants.atten_factorx = atten_factorx;
constants.atten_factory = atten_factory;
% 분자빔 초기조건
constants.jm_data           = jm_data;
constants.x_position_fwhm   = x_position_fwhm;
constants.Dye_sigma          = Dye_sigma;
constants.Dye_tau            = Dye_tau;
constants.vx_caused_tilted   = vx_caused_tilted;
constants.x_init_vel_sigma   = x_init_vel_sigma;
constants.y_init_vel_sigma   = y_init_vel_sigma;
constants.z_init_vel_off     = z_init_vel_off;
constants.z_init_vel_sigma   = z_init_vel_sigma;
constants.Ij_distribution    = Ij_distribution;
constants.Ips_z              = Ips_z;
constants.Ips_y              = Ips_y;
% 레이저 초기조건
constants.lambda = lambda;
constants.w0z    = w0z;
constants.w0y    = w0y;
constants.tau    = tau;
constants.atten1 = atten1;
constants.atten2 = atten2;
constants.tstep  = tstep;
constants.t_ini  = t_ini;
constants.tend   = tend;

out_file = fullfile(root_dir, sprintf('result_%s.mat', strrep(jm_data,'.txt','')));
save(out_file, 'Cell_velocity', 'I1_list', 'I2_list', 'constants');
fprintf('저장 완료: %s\n', out_file);

fprintf('\n=== 완료 '); toc
