function run_all(override)
% run_all(override) — PSW Dispersion Simulation
%
% override (선택): STN_constants.mat 기본값에서 바꿀 필드만 담은 struct
%   예) run_all(struct('jm_data','alpha_T200.txt','MAX',1000000))
%   인자 없이 run_all() 호출하면 기본값으로 실행
%
% 병렬 전략: parfor를 sub-batch 단위로 적용
%   MAX=200,000 × 8조건 → sub-batch 20,000 × 80태스크
%   8 workers (메모리 대역폭 포화점, 실측 최적)
%   실측 시간: ~181초 (기존 순차 ~960초 대비 5.3배 빠름)

if nargin < 1 || isempty(override)
    override = struct();
end

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

% 워커당 스레드 수 설정: 8 workers × 3 threads = 24 코어 활용
spmd
    maxNumCompThreads(3);
end


%% ── 1. 상수 로드 ─────────────────────────────────────────────────────────────
load(fullfile(root_dir, 'STN_constants_fullpaper.mat'));
MAX = 200000;   % 조건당 총 입자 수 (STN_constants의 값 overwrite)

% override struct의 필드를 로컬 워크스페이스에 덮어쓰기
% label은 출력 파일명용이므로 별도 추출 후 제외
if isfield(override, 'label')
    result_label = override.label;
else
    result_label = '';
end
override_fields = fieldnames(override);
for fi = 1:numel(override_fields)
    if ~strcmp(override_fields{fi}, 'label')
        eval(sprintf('%s = override.(''%s'');', override_fields{fi}, override_fields{fi}));
    end
end
clear override_fields fi

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

if numel(I1_list) ~= numel(I2_list)
    error('I1_list와 I2_list 크기가 다릅니다: I1=%d개, I2=%d개', numel(I1_list), numel(I2_list));
end

I1_vals = max_I * I1_list(:);   % n_groups×1
I2_vals = max_I * I2_list(:);   % n_groups×1

n_groups   = numel(I2_list);
batch_size = 100000;                          % sub-batch당 입자 수 (GPU VRAM 또는 캐시 기준)
n_sub      = ceil(MAX / batch_size);         % 조건당 sub-batch 수 (예: 10)
n_tasks    = n_groups * n_sub;               % 전체 parfor 태스크 수 (예: 80)

fprintf('[2/2] Time integration 시작...\n');
fprintf('      조건: %d × sub-batch: %d × 입자: %d = 총 태스크 %d\n', ...
    n_groups, n_sub, batch_size, n_tasks);
fprintf('      %d workers → ceil(%d/%d)=%d 라운드 예상\n', ...
    n_workers, n_tasks, n_workers, ceil(n_tasks/n_workers));

%% ── 6. parfor: sub-batch 단위 병렬 시뮬레이션 ───────────────────────────────
% broadcast 파라미터를 struct로 패킹하여 헬퍼 함수에 전달
P.x_position_fwhm  = x_position_fwhm;
P.Dye_sigma         = Dye_sigma;
P.Dye_tau           = Dye_tau;
P.vx_caused_tilted  = vx_caused_tilted;
P.x_init_vel_sigma  = x_init_vel_sigma;
P.y_init_vel_sigma  = y_init_vel_sigma;
P.z_init_vel_off    = z_init_vel_off;
P.z_init_vel_sigma  = z_init_vel_sigma;
P.Ij_distribution   = Ij_distribution;
P.Ips_z             = Ips_z;
P.Ips_y             = Ips_y;
P.lambda            = lambda;
P.w0z               = w0z;
P.w0y               = w0y;
P.tau               = tau;
P.tstep             = tstep;
P.t_ini             = t_ini;
P.tend              = tend;

results = cell(n_tasks, 1);

parfor task = 1:n_tasks
    results{task} = run_parfor_task(task, n_sub, batch_size, jm_data_path, ...
        linearly_2d, nj_lin, nm_lin, I1_vals, I2_vals, P);
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

if ~isempty(result_label)
    out_name = result_label;
else
    out_name = strrep(jm_data, '.txt', '');
end
out_file = fullfile(root_dir, 'result', sprintf('result_%s.mat', out_name));
save(out_file, 'Cell_velocity', 'I1_list', 'I2_list', 'constants');
fprintf('저장 완료: %s\n', out_file);

fprintf('\n=== 완료 '); toc

end
