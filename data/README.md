# 데이터 안내

이 프로젝트가 사용하는 외부 데이터 소스, 각 폴더의 역할, 저장소 포함 여부를 정리한다.

## 폴더 구조와 각 폴더의 의미

| 폴더 | 의미 |
|---|---|
| `data/raw/` | 원본 데이터. 외부에서 받은 그대로, 가공하지 않음 |
| `data/external/` | 용량 문제로 저장소에 포함하지 않는 대용량 외부 참조 데이터 (ClinVar) |
| `data/processed/` | `src/`의 앞 단계 스크립트가 만들고, 뒤 단계 스크립트가 다시 입력으로 읽는 "기계용" 중간 산출물 |
| `results/tables/` | 더 이상 다른 스크립트의 입력으로 쓰이지 않는, 사람이 검토하도록 만든 최종/해석용 결과 |

`data/processed/`와 `results/tables/`를 나누는 기준은 하나다 — **뒤 단계 파이프라인이 그 파일을
다시 읽는가**. 읽으면 `processed/`, 안 읽으면(그 자체로 최종 결과면) `results/tables/`.
어떤 스크립트가 무엇을 읽고 쓰는지는 [`docs/workflow.md`](../docs/workflow.md)에 단계별로 정리돼 있다.

**예외**: `data/processed/entrez_symbol_map.json`(Entrez 유전자 ID → 심볼 매핑)은 최초 커밋부터
존재했지만 현재 `src/`의 어떤 스크립트도 이 파일을 직접 만들거나 읽지 않는다 — 같은 매핑을
`src/kegg_analysis_utils.py`의 `build_gene_symbol_map()`이 원본 flat file에서 메모리 상으로
다시 계산한다. 삭제하지 않고 참고용 레거시 산출물로 보존했다.

---

## KEGG

**사용 파일** (`data/raw/kegg/`): `dgroup.txt`, `disease.txt`, `drug .txt`, `network.txt`, `variant.txt`
— KEGG 고정폭 flat file 포맷(`ENTRY ... ///`)의 원본 텍스트.

**이 프로젝트에서의 역할**: 질병(Disease) · 약물(Drug) · 약물군(DGroup) · 신호전달/질병
네트워크(Network) · 유전자 변이(Variant) 5개 엔트리 타입과, 엔트리 간 상호 참조 관계를
구조화된 CSV/edge list로 변환하는 원재료. `src/01_build_kegg_tables.py`가 이 5개 파일을
파싱해 `data/processed/kegg_*.csv`, `kegg_relations.csv`를 만든다.

**저장소 포함 여부**: 포함됨 (5개 파일 합쳐 약 35MB, git 추적).

**라이선스/재배포 관련 유의사항**: KEGG는 [자체 이용 약관](https://www.kegg.jp/kegg/legal.html)이
있는 데이터베이스다. 웹사이트를 통한 개인/학술적 열람은 자유롭지만, flat file을 대량으로
재배포하거나 상업적으로 이용하는 것은 별도 라이선스(Pathway Solutions/KEGG FTP 구독)가
필요할 수 있다. 이 저장소는 개인 포트폴리오·교육 목적의 비상업적 분석 예시이며, 저장소를
공개로 유지할 계획이라면 최신 KEGG 약관을 다시 확인하는 것을 권장한다.

---

## 가족 SNP 데이터셋 (Family SNP dataset)

**출처**: 공개 Kaggle 가족 유전체 데이터셋(23andMe 원시 genotype 내보내기 형식). 실제 살아있는
개인을 특정할 수 있는 이름 등은 포함돼 있지 않고, Father/Mother/Child 1/Child 2/Child 3의
5인 가족 구성으로 배포된 공개 데이터다.

**파일** (`data/raw/family_snp/`):
- `Father Genome.csv`, `Mother Genome.csv`, `Child 1/2/3 Genome.csv` — 실제 분석에 쓰는 가족
  구성원 5명의 원본 SNP 파일
- `Family Genome.csv` — Child 3와 바이트 단위로 동일한 중복 파일(원본 데이터셋에 그대로
  포함돼 있어 삭제하지 않고 원형 보존, 분석에는 사용하지 않음)
- `genome_zeeshan_usmani.csv` — 가족 구성원이 아닌, 데이터셋 제공자가 함께 배포한 별도의
  공개 참조 게놈(분석 대상 아님)
- `genome_file_description.csv` — 23andMe가 원본 내보내기 파일 머리말에 남긴 안내문(빌드
  버전, 컬럼 설명 등)

**원본 SNP/genotype 포맷**: 헤더 없는 4컬럼 CSV — `rsid, chromosome, position, genotype`.
`genotype`은 2글자 대립유전자 쌍(예: `AG`) 또는 무호출(`--`)이다. 기준 빌드는 **GRCh37
(Annotation Release 104)**, 23andMe가 2017년경 생성한 원시 데이터다.

**용도 제한**: 23andMe 원본 안내문에 명시된 대로 이 데이터는 "연구·교육·정보 제공 목적으로만
적합하며, 의료 또는 기타 용도로 사용할 수 없다." 이 프로젝트의 모든 분석은 그 제한을
그대로 따르는 가설 생성용 스크리닝이며, 임상적 해석으로 확장하지 않는다. 이 저장소의
[MIT License](../LICENSE)는 `src/`의 코드에만 적용되며, 이 데이터 자체의 라이선스를
대체하지 않는다.

**저장소 포함 여부**: 포함됨 (5개 원본 파일 합쳐 약 75MB, git 추적). `*Copy.csv` 형태의
바이트 단위 중복 사본은 `.gitignore`로 제외했다(원본 데이터셋에 딸려온 백업 사본으로 추정,
git 이력에 없던 파일이라 삭제해도 데이터 손실이 없음).

---

## NCBI ClinVar

**사용 파일**: `variant_summary.txt` (탭 구분 텍스트, gzip 압축본 `variant_summary.txt.gz`로 배포)

**이 프로젝트에서의 역할**: 가족 rsid를 실제 임상 유의성(ClinicalSignificance)과 rsid 단위로
연결하는 유일한 외부 임상 데이터베이스. `src/05_clinvar_annotation.py`가 GRCh37 빌드 행만
필터링하고, 가족 5명의 rsid 합집합(1,117,586개)과 매칭되는 행만 스트리밍으로 추출한다.

**genome build**: GRCh37 (가족 23andMe 원본 데이터의 기준 빌드와 동일하게 맞춤).

**저장소 포함 여부**: **포함되지 않음.** 압축 해제 시 약 3.7GB, 압축 상태로도 약 420MB로,
git 저장소에 넣기에는 너무 크다 (`data/external/clinvar/`는 `.gitignore`로 제외).

**재현 방법**:
1. [NCBI ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)에서
   `variant_summary.txt.gz`를 내려받는다.
2. `data/external/clinvar/variant_summary.txt.gz` 경로에 저장한다 (파일명 그대로).
3. `python src/05_clinvar_annotation.py`를 실행한다. 파일이 없으면 스크립트가 위와 동일한
   안내 메시지를 출력하고 즉시 종료한다.

ClinVar는 릴리스마다 내용이 갱신되므로, 이 저장소에 커밋된 `family_clinvar_matches.csv` 등의
결과는 프로젝트 작업 시점(2026-08)의 ClinVar 릴리스를 기준으로 한다 — 최신 파일로 재현하면
매칭 건수가 소폭 달라질 수 있다.

---

## 요약: 원본 vs 외부 vs 가공 vs 결과

| 구분 | 위치 | 예시 |
|---|---|---|
| 원본 데이터(raw) | `data/raw/kegg/`, `data/raw/family_snp/` | `dgroup.txt`, `Father Genome.csv` |
| 외부 참조 데이터(external, 저장소 미포함) | `data/external/clinvar/` | `variant_summary.txt.gz` |
| 가공 중간 산출물(processed) | `data/processed/` | `kegg_relations.csv`, `family_clinvar_matches.csv` |
| 분석 결과(results) | `results/tables/`, `results/figures/`, `results/case_studies/` | `family_pharmacogenomics.csv`, CFH 사례 연구 |
