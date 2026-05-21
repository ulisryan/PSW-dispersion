% run_batch.m — 복수 조건 배치 실행
%
% 사용법:
%   각 param_set은 STN_constants.mat 기본값에서 바꿀 필드만 지정한 struct.
%   'label' 필드(선택): 결과 파일명 suffix. 없으면 jm_data 이름 사용.
%
% 예시:
%   struct('jm_data','alpha_T300.txt')
%     → result_alpha_T300.mat 저장
%   struct('jm_data','alpha_T200.txt','z_init_vel_off',550,'label','T200_v550')
%     → result_T200_v550.mat 저장, z_init_vel_off만 변경

clear; clc;

%% ── 파라미터 세트 정의 ────────────────────────────────────────────────────────
% 아래 param_sets 셀 배열에 원하는 조건을 추가하세요.
% 지정하지 않은 필드는 모두 STN_constants.mat 기본값 사용.

vx_values = [8:2:10];
param_sets = cell(1, numel(vx_values));
for k = 1:numel(vx_values)
% for k = 5:6
    param_sets{k} = struct( ...
        'vx_caused_tilted', vx_values(k), ...
        'label', sprintf('20260515_zoom_%02d', k+4), ...
        'I1_list', [1.07:0.0005:2.13], ...
        'I2_list', [1.07:0.0005:2.13] ...
    );
end

%% ── 배치 실행 ────────────────────────────────────────────────────────────────
n_sets = numel(param_sets);
fprintf('총 %d 세트 실행 시작\n\n', n_sets);
batch_tic = tic;

for k = 1:n_sets
    fprintf('========== 세트 %d / %d ==========\n', k, n_sets);
    disp(param_sets{k});   % 현재 세트 내용 출력
    run_all(param_sets{k});
    fprintf('\n');
end

fprintf('===== 모든 배치 완료 (총 ');
toc(batch_tic);
