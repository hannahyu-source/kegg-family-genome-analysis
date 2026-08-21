# 작업 로그 — KEGG × 가족 유전체 분석 프로젝트

> **역사적 기록**: 이 문서는 최초 분석 세션(2026-08-21)의 작업 기록을 그대로 보존한 것이다.
> 이 문서에 나오는 `KEGG 데이터/script/`, `가족 유전체(SNP) 데이터/output/` 같은 경로는
> 그 세션 당시의 폴더 구조를 가리키며, 이후 저장소를 포트폴리오 구조(`src/`, `data/`,
> `results/`)로 재구성하면서도 원래 기록을 왜곡하지 않기 위해 수정하지 않았다. 현재 구조
> 기준의 파일 위치는 [`workflow.md`](workflow.md)를 참고할 것.

기간: 2026-08-21 (단일 세션)
범위: `KEGG 데이터/`, `가족 유전체(SNP) 데이터/` 두 폴더 (재구성 이전 구조)

최종 산출물 3종:
- 보고서: **Pathway to Pedigree** — https://claude.ai/code/artifact/fbb4aa8c-875b-4c25-a8ac-3972394a80e1
- 시각화: **KEGG Relations Atlas** — https://claude.ai/code/artifact/5dab17ec-13a4-4f0a-9203-ca4466ab5426
- 시각화: **Drug Taxonomy** — https://claude.ai/code/artifact/9305b69f-25a7-4c4f-a1b9-fa129b9bbb85

---

## 0. 폴더 정리

- `KEGG 데이터/`, `가족 유전체(SNP) 데이터/` 두 폴더 생성
- 루트에 섞여 있던 파일을 확장자로 분리: `*.txt` → `KEGG 데이터/`, `*.csv` → `가족 유전체(SNP) 데이터/`
- 두 폴더 모두 `data / docs / output / script` 4분할 구조로 정리, 원본 파일은 각 `data/`로 이동

---

## 1. KEGG 데이터 구조화

**스크립트**: `KEGG 데이터/script/kegg_flatfile_parser.py`(공용 파서), `build_kegg_tables.py`

- KEGG flat file(`ENTRY … ///`, 12칸 고정폭 필드) 공용 파서 작성
- 5개 원본 파일 → 구조화 CSV 5개 + 관계 엣지 테이블 1개

| 원본 | 엔트리 수 | 출력 |
|---|---|---|
| dgroup.txt | 2,429 | kegg_dgroup.csv |
| disease.txt | 2,633 | kegg_disease.csv |
| drug .txt | 12,115 | kegg_drug.csv |
| network.txt | 1,417 | kegg_network.csv |
| variant.txt | 982 | kegg_variant.csv |

- 엔트리 간 관계 **69,546건**을 `kegg_relations.csv`로 통합 (bare ID·`[TAG:ID]` 형태 참조 정규식 추출)
- `variant.txt` GENE 필드만 따로 뽑아 `variant_genes.csv` / `variant_genes_unique.csv`(909개 질병연관 유전자) 생성 (`extract_variant_genes.py`)

**버그 발견·수정**: REMARK 필드를 정규식 lookahead로 자르다가 다음 항목까지 같이 잘려나가는 문제 발견. 줄 단위 key:value 파싱으로 교체 후 `drug→compound` 283→1,693건, `drug→dgroup` 1,771→4,345건으로 정상화.

---

## 2. 관계 그래프 분석·시각화

**스크립트**: `analyze_kegg_relations.py`

- `kegg_relations.csv`를 pandas groupby로 집계 (networkx 불필요)
- 출력: `kegg_graph_type_summary.csv`, `kegg_graph_top_hubs.csv`, `kegg_graph_ego_network.json`
- 허브: dgroup "CYP3A/CYP3A4 substrate"(1,024개 약물), disease는 간암(H00048, degree 139)이 최다
- `network.txt`/`disease.txt` GENE 필드에서 Entrez ID → 심볼 매핑 테이블 별도 구축(5,953개)
- 간암 1-hop 관계망(140노드: gene 22·drug 12·variant 18·network 87) 시각화 → **KEGG Relations Atlas** 아티팩트
  - 동심원(orbit) 레이아웃 채택 (물리 시뮬레이션 시도했으나 스타 토폴로지라 부적합 판단)
  - 색상: dataviz 팔레트 중 all-pairs CVD 검증 통과한 3색(blue/aqua/orange)만 사용, 나머지는 회색 도형(원/마름모)으로 구분

---

## 3. 가족 SNP 교차분석

**스크립트**: `가족 유전체(SNP) 데이터/script/family_snp_cross_analysis.py`, `family_clinvar_match.py`, `family_clinvar_kegg_crossvalidate.py`

### 3-1. 데이터 정리
`data/` 13개 파일 중 실제 가족 구성원은 5명(Father/Mother/Child 1/2/3)뿐. `*-Copy.csv`(중복), `Family Genome.csv`(Child 3와 완전 동일), `genome_zeeshan_usmani.csv`(가족 아닌 공개 참조 게놈) 제외.

### 3-2. 기초 통계
`family_snp_summary.csv` — 구성원별 총 SNP(60만~63만), 결측률(1.1~2.5%), 이형접합률. Child 2/3이 나머지보다 뚜렷이 낮음(16.7~17.3% vs 28.5~28.9%) — 유전체칩 버전 차이로 추정되는 QC 관찰.

### 3-3. 1차 시도 (좌표 없이 가능한 범위)
KEGG variant.txt가 rsid가 아니라 OMIM 번호 기준이라 직접 매칭 불가 → `kegg_disease_gene_panel.csv`(909개 유전자, 질병 연결 수 기준 순위)만 우선 작성.

### 3-4. 사용자 제공 파일 검토
`clinvar_ml_balanced.csv`(640,631행) 검토 결과 rsid/유전자/좌표 없는 병원성 판정 ML 피처 테이블로 확인, 매칭 불가 판정.

### 3-5. 실제 ClinVar 매칭
NCBI ClinVar FTP `variant_summary.txt.gz`(442MB, 9,044,811행) 다운로드 → GRCh37 필터링(4,541,938행) → 가족 rsid(1,117,586개 union) 스트리밍 매칭 → **59,501건** 매칭, 그중 55,730건 SNV 대립유전자 비교 가능.
출력: `family_clinvar_matches.csv`, `family_clinvar_summary.csv`

**발견**: Pathogenic/Likely pathogenic(비상충) 실제 보유 4건 — CFH×2(rs460897 ★★10명, rs1061170 ★1명=흔한 Y402H 다형성), LMNA×1(rs58922911, ★0), UGT1A1×1(rs3755319, ★0)

### 3-6. KEGG 교차검증
ClinVar 병원성 4건을 `kegg_disease_gene_panel.csv`(명시적 DISEASE: H##### 링크 보유 872개 유전자)와 join → **CFH 2건만 통과**, LMNA·UGT1A1은 KEGG에 명시적 질병 링크가 없어 탈락.
출력: `family_clinvar_kegg_crossvalidated.csv`

---

## 4. 질병 중심 분석

**스크립트**: `KEGG 데이터/script/disease_centered_analysis.py`

- 질병의 (직접 유전자 + 관련 network 유전자) = 병리경로 유전자 집합 → 그 유전자를 표적하는 미승인 약물 = 재창출 후보
- network↔disease 링크는 network.txt DISEASE 필드(bare H#####) 방향만 신뢰 가능 (disease.txt NETWORK 필드는 nt0xxxx CLASS 코드라 별개 네임스페이스 — 직접 연결 불가, 데이터 구조 확인 후 판단)
- 전체 2,633개 질병 중 673개에서 재창출 후보 발견, 총 5,890건
  - 검증 예: 간암 → PI3K/mTOR 억제제(MTOR/PIK3CA/PIK3CB/PIK3CD 표적) — 실제 임상연구 표적과 일치
- 동반질환(공유 network/gene) 상위 300쌍: Burkitt/Hodgkin lymphoma/비인두암(EBV 연관암, 30 networks), 비대성/확장성 심근병증(사르코미어 유전자 공유, 18 genes), 뇌하수체선종↔쿠싱증후군(실제 인과관계) 등 — 전부 실제 의학 지식과 부합
출력: `disease_pathway_genes.csv`, `drug_repurposing_candidates.csv`, `disease_comorbidity_pairs.csv`

---

## 5. 가족 약물유전체(PGx) 분석

**스크립트**: `가족 유전체(SNP) 데이터/script/family_pharmacogenomics.py`

- `family_clinvar_matches.csv`에서 ClinicalSignificance="drug response" 실제 보유 60건(유전자 26개) 추출
- KEGG dgroup 이름(예: "CYP3A4 substrate") + PhenotypeList의 "약물명 response" 패턴 → KEGG drug 엔트리 직접 매칭(53/60건 성공)
- **전 가족(5명) 공통 변이**: VKORC1(warfarin), CYP3A5(tacrolimus), UGT1A1(irinotecan) — 전부 실제 CPIC/FDA 임상 가이드라인 마커
출력: `family_pharmacogenomics.csv`, `family_pharmacogenomics_summary.csv`

---

## 6. 유전자 영역 스크리닝 (독립 검증)

**스크립트**: `family_kegg_gene_screening.py`

- KEGG 909개 유전자 × family_clinvar_matches.csv를 유전자 심볼로 join
- **버그 발견·수정**: 처음엔 "테스트된 위치가 있음"과 "실제 대립유전자 보유"를 구분 못해 BRCA1/BRCA2 등 병원성 기록이 많은 유전자가 잘못 상위에 뜸 → `carries_alt_allele` 필터 추가로 수정
- 결과: 846개 유전자에 테스트 위치 존재, 601개는 실제 보유 확인(2,069 rsid) — Pathogenic 플래그는 CFH/UGT1A1/LMNA로 Phase 3-5와 정확히 일치(독립 재현 확인)
출력: `family_kegg_gene_screening.csv`

---

## 7. CFH 집중 프로파일

**스크립트**: `cfh_family_disease_profile.py`

- Phase 3에서 확인한 CFH 실제 보유 변이를 Phase 4 엔진(재창출 후보 + 동반질환)에 대입
- CFH 연결 5개 질병(대체 보체 경로 결함/황반변성/비정형 용혈성요독증후군/기저층 드루젠/C3 사구체병증) 전부 보체계 질환으로 강하게 클러스터링
- 재창출 후보: Pegcetacoplan, Iptacopan, Danicopan — 실제 FDA 승인/후기임상 보체억제제와 일치
출력: `cfh_family_disease_profile.csv`

---

## 8. 약물 데이터 활용

**스크립트**: `KEGG 데이터/script/atc_classification_stats.py`, `dgroup_class_tree.py`

- ATC 코드 1단계(해부학적 주분류) 기준 치료영역별 통계: N 신경계(885) > A 소화기·대사(784) > J 항감염제(773) > C 심혈관계(705)
- `dgroup.txt` CLASS 필드 파싱 → 다중부모 DAG 구조 확인(텍스트 루트 카테고리 → dgroup 체인, 최대 3~4단계)
- 29개 최상위 카테고리: "Metabolizing enzyme substrate"(1,155개 약물) 최다, "Neuropsychiatric agent"(1,020개)가 ATC N분류와 교차검증됨
- 펼침 트리(`<details>`) 인터랙티브 시각화 → **Drug Taxonomy** 아티팩트
출력: `atc_top_level_summary.csv`, `atc_subgroup_summary.csv`, `dgroup_class_edges.csv`, `dgroup_class_rollup.csv`

---

## 산출물 전체 목록

### `KEGG 데이터/output/`
```
kegg_dgroup.csv  kegg_disease.csv  kegg_drug.csv  kegg_network.csv  kegg_variant.csv
kegg_relations.csv  (69,546 edges)
variant_genes.csv  variant_genes_unique.csv
kegg_graph_type_summary.csv  kegg_graph_top_hubs.csv  kegg_graph_ego_network.json
disease_pathway_genes.csv  drug_repurposing_candidates.csv  disease_comorbidity_pairs.csv
atc_top_level_summary.csv  atc_subgroup_summary.csv
dgroup_class_edges.csv  dgroup_class_rollup.csv
```

### `가족 유전체(SNP) 데이터/output/`
```
family_snp_summary.csv  family_snp_by_chromosome.csv
kegg_disease_gene_panel.csv  (909 genes)
family_clinvar_matches.csv  (59,501 rows)  family_clinvar_summary.csv
family_clinvar_kegg_crossvalidated.csv  (2 rows)
family_pharmacogenomics.csv  (60 rows)  family_pharmacogenomics_summary.csv
family_kegg_gene_screening.csv  (601 genes)
cfh_family_disease_profile.csv
```

### 스크립트
```
KEGG 데이터/script/
  kegg_flatfile_parser.py  build_kegg_tables.py  extract_variant_genes.py
  analyze_kegg_relations.py  disease_centered_analysis.py
  atc_classification_stats.py  dgroup_class_tree.py

가족 유전체(SNP) 데이터/script/
  family_snp_cross_analysis.py  family_clinvar_match.py
  family_clinvar_kegg_crossvalidate.py  family_pharmacogenomics.py
  family_kegg_gene_screening.py  cfh_family_disease_profile.py
```

---

## 참고 / 제약사항

- 참조 데이터: NCBI ClinVar `variant_summary.txt` (GRCh37, 2026-08 릴리스)
- 가족 SNP 원본: 23andMe raw genotype, GRCh37(Annotation Release 104) 기준
- 이 프로젝트의 모든 "재창출 후보"·"동반질환"·"PGx 매칭"·"CFH 프로파일"은 **가설 생성·참고용 스크리닝 결과**이며 의학적 진단·처방 근거가 아님. 실제 임상 적용 전에는 반드시 의사·약사·유전상담사와 상의할 것.
