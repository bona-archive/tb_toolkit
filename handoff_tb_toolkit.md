# Handoff: tb_toolkit — 응집물질물리 Tight-Binding 라이브러리

## 컨텍스트

사용자는 응집물질물리 tight-binding 계산에 쓰이는 MATLAB 유틸리티 함수들을 **범용 외부 라이브러리**로 정리하는 중이다. 원래 `C:\Users\hsb83\Documents\MATLAB\kondo_kagome`에 kagome 전용으로 작성된 함수들을 출발점으로 삼았다.

**중요: Claude의 역할은 튜터이지 코드 작성자가 아니다.** CLAUDE.md에 명시되어 있음. 사용자가 직접 코드를 작성하고, Claude는 질문/피드백/힌트로 유도해야 한다. 완성된 코드를 직접 써주면 안 된다. 항상 한국어로 응답.

---

## 현재 상태

### 라이브러리 위치
```
C:\Users\hsb83\Documents\MATLAB\tb_toolkit\
  setup.m
  +tb\
    lat_2d.m
    prim2reci.m
    reci2BZ.m
    distmat.m
    Draw_coupling_graph.m
    BZmask.m
    hexBZverts.m
    kpath.m         ← 사용자가 이번 세션에 직접 작성
    blochham.m      ← 사용자가 이번 세션에 직접 작성
    dist_decay_hops.m ← 사용자가 이번 세션에 직접 작성
```

### 이번 세션에서 완성된 함수들

**`tb.kpath`**
- 입력: `kpts` (cell array), `labels` (TeX strings cell array), `N_pts`
- 출력: `k_path`, `k_dist`, `ticks`, `labels`
- 경계점 중복 제거 로직 포함

**`tb.blochham`**
- 입력: `k_vec (1,2)`, `basis (:,2)`, `latvecs (2,2)`, `hops (:,4)`, `onsite (:,1)`
- `hops` 형태: `[i, j, dist, t]`
- Bloch 위상 gauge: `d_vec = R + b_j - b_i` (orbital exact position gauge)
- 출력: `Hk` (n_orb × n_orb Hermitian)

**`tb.dist_decay_hops`**
- distance-decay 호핑 편의 함수: `t * exp(-gamma * dist)`
- 입력: `basis`, `latvecs`, `neighbor_order`, `t`, `gamma`, `options.doPlot`
- 내부적으로 3×3 슈퍼셀 생성 후 distmat으로 이웃 탐색

### 테스트 스크립트 (작동 확인됨)
honeycomb 2-band 모델에서 `tb.blochham` + `tb.kpath`로 밴드 구조 계산 성공.

---

## 다음 세션의 주요 작업: 네임스페이스 재구조화

사용자가 결정한 Option A 구조:

```
tb_toolkit/
  +lat/         ← 격자, 거리
    lat_2d.m
    distmat.m

  +bloch/       ← k-space 계산
    kpath.m
    blochham.m

  +obc/         ← OBC/리본 계산
    dist_decay_hops.m

  +draw/        ← 시각화
    Draw_coupling_graph.m

  +to/          ← 변환 (transformation)
    prim2reci.m
    reci2BZ.m

  +mask/        ← BZ 마스킹
    BZmask.m
    hexBZverts.m

  +Presets/     ← 잘 알려진 격자 기본값 (이름 미확정)
    graphene.m
    kagome.m
    triangular.m
    square.m
    lieb.m

  +SymmOp/      ← 대칭 연산
    Rotation.m
    Mirror.m
    Parity.m
    Translation.m

  setup.m
```

**미결 사항:**
- `+Presets` 이름이 적절한지 사용자가 고민 중
- `+to` 네임스페이스 이름이 직관적인지 재검토 필요 (`+reci`? `+reciprocal`?)
- 기존 `+tb` 함수들을 새 네임스페이스로 이동/리네임

---

## 다음 세션에서 할 일 (우선순위 순)

1. **네임스페이스 재구조화 실행** — 파일들을 새 폴더 구조로 이동
2. **`obc.edgeweight` 구현** — `tight_binding_OBC_y.m` 99~119번 줄 기반. 파동함수에서 edge localization weight 계산
3. **`+Presets` 구현** — graphene, kagome 등 표준 격자 파라미터 반환 함수
4. **테스트** — honeycomb 모델에서 새 네임스페이스로 전체 워크플로 검증

---

## 참고 파일

- 원본 코드: `C:\Users\hsb83\Documents\MATLAB\kondo_kagome\`
- Wannier topology 프로젝트: `C:\Users\hsb83\Documents\Matlab\project_Wannier_topology\matlab_code\`
  - `tight_binding_OBC_y.m` — OBC 리본, edge state 계산 (obc 모듈 참고용)
  - `tight_binding_primitive_cell_band.m` — Bloch H(k) 계산 참고
  - `zak_phase.m` — Wilson loop / Zak phase (미래 +bloch 확장 참고)

---

## 주의사항

- `tb.blochham`은 `hops`에 같은 `(i,j,dist)` 조합이 중복되면 `match`가 여러 행 반환 → 문제 발생. `dist_decay_hops`는 `unique('rows')`로 처리함.

---

## Suggested Skills

- `/caveman` — 사용자가 토큰 절약 원할 때
- `/handoff` — 다음 세션 마무리 시
