# Validation Strategy

이 프로젝트가 결과를 어떻게 검증했는지, 그리고 검증하지 못한 부분은 어디인지 명시한다.
아래 검증 중 다수는 이번 저장소 재구성 작업 중 파이프라인을 실제로 재실행해서 확인했다
(2026-08-21) — 재실행 로그와 수치는 [`development_log.md`](development_log.md) 및 각
스크립트 실행 출력과 대조했다.

## 데이터 무결성 검사

- **파일 존재 확인**: 각 스크립트는 `pathlib.Path` 기반 상대 경로로 입력 파일을 읽으며,
  ClinVar처럼 저장소에 없는 대용량 외부 파일은 없을 경우 다운로드 안내와 함께
  `FileNotFoundError`를 던지도록 만들어 "조용히 잘못된 결과를 내는" 상황을 방지한다
  (`src/05_clinvar_annotation.py`).
- **입력 행 수**: KEGG 5개 flat file → 엔트리 수 dgroup 2,429 / disease 2,633 / drug 12,115 /
  network 1,417 / variant 982건, 관계 69,546건으로 파싱됨을 재실행으로 재확인.
- **중복 rsid**: KEGG 매칭 대상 rsid는 가족 5명의 rsid를 `dict`(rsid→genotype)로 만들어
  합집합을 구하므로 구성원 내부 중복은 자동으로 제거된다(같은 rsid가 파일에 두 번 나오면
  마지막 값으로 덮어써짐 — 23andMe 원본 형식상 발생 가능성은 낮지만 별도 중복 카운트 로직은
  없음, 향후 개선 여지로 [`limitations.md`](limitations.md)에 기록).
- **결측 genotype**: `family_snp_summary.csv`의 `no_call_count`/`no_call_rate`로 구성원별
  집계(1.1%~2.5%). 별도 대체(imputation) 없이 그대로 보고한다.
- **염색체 값**: `family_snp_by_chromosome.csv`로 구성원 x 염색체 조합별 SNP 개수를 집계 —
  각 구성원 25개의 distinct chromosome 값(상염색체 1–22, X, Y, MT)을 확인.
- **가족 표본 수**: 5명(Father/Mother/Child 1/2/3) — 13개 원본 파일 중 중복/비가족 파일을
  제외하고 명시적으로 고정.
- **병합 카디널리티**: ClinVar 매칭은 `rsid` 단일 키로 가족 파일과 조인하며, ClinVar 쪽에
  동일 rsid가 여러 VariationID로 여러 번 등록된 경우 1:N으로 매칭되어 매칭 결과 행 수가
  가족 고유 rsid 수보다 많아질 수 있다(실제로 매칭 rsid 59,501건 vs 가족 고유 rsid
  1,117,586개 — ClinVar 쪽에 존재하는 것만 남으므로 매칭 후 행 수가 더 적다).

## ClinVar 매칭 검증

- **유효 rsid**: `RS# (dbSNP) != "-1"`인 행만 사용(RS#가 없는 항목 제외).
- **genome build**: `Assembly == "GRCh37"` 필터를 명시적으로 적용, 가족 데이터의 기준 빌드와
  일치시킴.
- **ClinVar 중복 레코드**: 같은 rsid가 여러 VariationID(서로 다른 변이 해석)로 등록될 수
  있음 — 이 프로젝트는 이를 제거하지 않고 모두 유지한다(한 rsid가 여러 임상적 해석을 가질
  수 있다는 사실 자체가 정보이기 때문). 다만 그만큼 "매칭 건수"는 "고유 rsid 수"가 아니라
  "rsid-VariationID 쌍의 수"에 가깝다는 점을 결과 해석 시 유의해야 한다.
- **ClinicalSignificance 해석**: 문자열 포함 검사로 `"Pathogenic"`을 병원성 계열로 분류하되
  `"Conflicting"`이 포함된 항목은 제외한다(제출자 간 의견 불일치 항목을 병원성으로
  단정하지 않기 위함).
- **review status**: `ReviewStatus`, `NumberSubmitters` 컬럼을 원본 그대로 보존해
  `family_clinvar_matches.csv`에 남긴다 — 예를 들어 CFH의 rs460897은 10명 제출자 합의
  (★★, 상대적으로 신뢰도 높음), rs1061170은 제출자 1명·무검토(★, 흔한 다형성 Y402H)로
  신뢰도 차이가 크다는 점을 원본 데이터 그대로 확인할 수 있다.
- **allele 모호성 한계**: SNV·단일 문자 Ref/Alt인 경우에만 `carries_alt_allele`를 계산하고,
  그 외(indel 등)는 `None`으로 남겨 "모른다"와 "보유하지 않는다"를 구분한다.

## KEGG 교차검증

KEGG는 ClinVar가 내린 병원성 판정을 "검증"하는 근거가 아니라, **서로 다른 큐레이션 과정을
거친 독립적인 지식 소스**로 사용한다. 교차검증이 뜻하는 바는 "이 유전자가 두 데이터베이스
모두에서 질병과 연결돼 있다"는 것뿐이며, "ClinVar가 병원성으로 판정한 그 변이의 정확한
기전을 KEGG가 재확인했다"는 뜻이 아니다. 실제로 가족이 보유한 ClinVar 병원성 변이 4건 중
KEGG 교차검증을 통과한 것은 CFH 2건뿐이었고(LMNA·UGT1A1 탈락), 이 탈락 자체가 "두 소스가
항상 일치하지는 않는다"는 것을 보여주는 유효한 결과로 취급한다.

## PGx 검증

**마커 식별(marker identification)**과 **임상 표현형/용량 해석(phenotype/dosing
interpretation)**을 명확히 구분한다:

- 이 프로젝트가 하는 것: ClinVar가 `"drug response"`로 분류한 변이를 가족이 보유하는지
  확인하고, KEGG dgroup/drug 이름과 텍스트 매칭한다.
- 이 프로젝트가 하지 않는 것: star-allele/haplotype 판정, 대사자 표현형(poor/intermediate/
  extensive/ultrarapid metabolizer) 분류, 실제 처방 용량 계산.
- **CPIC/FDA와의 관계**: VKORC1(warfarin)·CYP3A5(tacrolimus)·UGT1A1(irinotecan)은 실제
  CPIC/FDA 가이드라인에 마커로 등재된 유전자라는 사실은 공개된 CPIC/PharmGKB 자료와 대조해
  수작업으로 확인했다 — 다만 이 프로젝트 코드 자체가 CPIC 가이드라인 API나 파일을
  프로그램적으로 조회하지는 않는다(수작업 대조이지 자동 통합이 아님이라는 점을 명시).

## KEGG 관계 그래프 검증

`kegg_graph_top_hubs.csv`의 허브 결과(예: 간암 H00048이 disease 중 최다 연결, dgroup
"CYP3A/CYP3A4 substrate"가 최다 연결 약물군)는 KEGG 원본 데이터를 pandas groupby로 집계한
것으로, 별도 그래프 라이브러리(NetworkX 등)를 쓰지 않고 관계형 집계만으로 재현 가능하다.
`disease_centered_analysis.py`가 만드는 "재창출 후보"(예: 간암 → PI3K/mTOR 억제제)는 이
프로젝트 실행 중 실제 임상시험에서 연구되는 표적과 일치함을 수작업으로 확인했다 — 이 역시
자동 검증이 아니라 사람이 결과를 보고 생물학적 타당성을 판단한 것이다(아래 "수작업 검토"
참고).

## 수작업 검토(manual review)

다음 항목은 결과가 나온 뒤 사람이 직접 생물학적 타당성을 검토했다:
- CFH 관련 5개 질병이 모두 보체계(complement system) 질환으로 클러스터링되는지
- CFH 재창출 후보(Pegcetacoplan, Iptacopan, Danicopan)가 실제 FDA 승인/후기임상 보체억제제와
  일치하는지
- 간암 재창출 후보(PI3K/mTOR 억제제 계열)가 실제 임상연구 표적과 일치하는지
- 동반질환 상위 쌍(Burkitt/Hodgkin lymphoma/비인두암의 EBV 연관, 비대성/확장성 심근병증의
  사르코미어 유전자 공유, 뇌하수체선종↔쿠싱증후군의 실제 인과관계)이 알려진 의학 지식과
  부합하는지
- PGx 공통 마커(VKORC1/CYP3A5/UGT1A1)가 실제 CPIC/FDA 등재 마커인지

이 검토는 "결과가 그럴듯해 보인다"는 사후 확인이며, 전문의·유전상담사에 의한 공식 임상
검토가 아니다.

## AI 생성/수정 코드 검증

이 저장소의 코드는 Claude Code의 도움을 받아 작성·리팩터링됐다(자세한 내용은
[`ai_assisted_workflow.md`](ai_assisted_workflow.md)). AI가 만들거나 수정한 코드는 다음
방법으로 검증했다:

- **입출력 대조**: 저장소 재구성 과정에서 모든 경로를 새 구조로 바꾼 뒤 파이프라인을
  실제로 재실행해, 재구성 전 커밋에 존재하던 수치(69,546건 관계, 1,117,586개 가족 rsid,
  5,890건 재창출 후보, CFH 교차검증 2건, PGx 공통 마커 3개 등)와 정확히 일치하는지 확인했다.
- **행 수 비교**: 각 스크립트 실행 로그의 건수를 `development_log.md`에 기록된 원래 수치와
  1:1로 대조했다.
- **스팟 체크**: CFH rs460897/rs1061170 보유자 목록, PGx 공통 유전자 목록처럼 구체적인 값을
  직접 눈으로 대조했다.
- **독립 데이터베이스 비교**: KEGG 교차검증(6단계)과 유전자 스크리닝(8단계)이 서로 다른
  경로로 CFH를 동일하게 가리키는지 확인 — 이는 "AI가 만든 두 스크립트가 서로 모순되지
  않는다"는 재현성 검증이기도 하다.
- **생물학적 타당성 검토**: 위 "수작업 검토" 항목 전체.
- **문서 내 다이어그램 검증**: README.md의 파이프라인 Mermaid 다이어그램은 GitHub가 클라이언트
  측 JavaScript로 렌더링하기 때문에 정적 HTML만 가져오는 도구(예: 단순 URL fetch)로는 실제
  렌더링 여부를 확인할 수 없다 — 그런 도구로 확인하면 다이어그램이 원문 텍스트 그대로
  보여 "깨졌다"는 오탐(false negative)이 나온다. 이를 가리기 위해 README.md에서 Mermaid
  소스를 그대로 추출해 로컬에 설치한 `@mermaid-js/mermaid-cli`(GitHub와 같은 계열의
  mermaid.js 엔진)로 직접 렌더링했고, 오류 없이 정상적으로 그림이 만들어짐을 확인했다.

리팩터링 과정에서 발견한 것: 파이프라인 스크립트 중 하나(`family_clinvar_match.py`, 현재
`src/05_clinvar_annotation.py`)가 이 세션 밖에서 만들어진 임시 피클 캐시 파일과 로컬 절대
경로에 의존하고 있어 원래 상태로는 다른 환경에서 재현이 불가능했다 — 이를 원본 SNP CSV에서
직접 계산하도록 수정했다.

이후 `data/external/clinvar/variant_summary.txt.gz`가 로컬에 준비된 김에, 900만 행 전체
ClinVar 스트리밍 재매칭을 실제로 처음부터 끝까지 재실행해 최종 검증했다: 9,044,810행 중
GRCh37 4,541,938행, 가족 rsid와 매칭 **59,501건**으로 기존 커밋과 정확히 일치했고
(`git status`로 바이트 단위 무변경 확인), allele 비교 가능 55,730건도 동일했다. 즉
`05_clinvar_annotation.py`는 부분 검증이 아니라 처음부터 끝까지 완전히 재현 가능함을
실제로 확인했다.
