# KEGG × 가족 유전체(SNP) 교차분석

KEGG(질병·유전자·경로·약물 관계형 데이터베이스)를 구조화하고, 가족 5명의 23andMe SNP
데이터를 NCBI ClinVar와 실제로 매칭해서 KEGG 지식과 교차검증한 프로젝트.

## 데이터

| 소스 | 내용 | 위치 |
|---|---|---|
| KEGG flat file | dgroup / disease / drug / network / variant (5개) | `KEGG 데이터/data/` |
| 가족 SNP | Father/Mother/Child 1·2·3 raw genotype (Kaggle 공개 데이터셋) | `가족 유전체(SNP) 데이터/data/` |
| NCBI ClinVar | `variant_summary.txt`(GRCh37, 재현 시 직접 다운로드 필요 — 용량 문제로 repo 미포함) | 외부 |

## 폴더 구조

```
KEGG 데이터/                    가족 유전체(SNP) 데이터/
├── data/    원본 flat file      ├── data/    원본 SNP CSV
├── docs/                        ├── docs/
├── output/  구조화 CSV/JSON     ├── output/  분석 결과 CSV
└── script/  파이썬 스크립트     └── script/  파이썬 스크립트
```

## 파이프라인

```
KEGG flat file 파싱 → 관계 그래프 분석 → 가족 SNP × ClinVar 매칭 → KEGG 교차검증
→ 질병 중심 분석(재창출/동반질환) → 가족 약물유전체(PGx) → 유전자 스크리닝(독립검증)
→ CFH 프로파일 → 약물 분류(ATC/dgroup 트리)
```

단계별 입출력과 재현 실행 순서는 **[workflow.md](workflow.md)** 참고.
전체 작업 기록(발견한 버그, 수치, 판단 근거)은 **[work_log.md](work_log.md)** 참고.

## 주요 발견

- KEGG 5개 파일 → **69,546건**의 엔트리 간 관계로 구조화
- 가족 rsid 1,117,586개 중 ClinVar와 **59,501건** 매칭
- 가족이 실제 보유한 병원성(Pathogenic) 변이 중 **CFH**(rs460897, rs1061170)만 KEGG
  자체 질병 데이터베이스로도 교차검증됨
- 전 가족 공통 약물유전체(PGx) 마커: **VKORC1**(warfarin), **CYP3A5**(tacrolimus),
  **UGT1A1**(irinotecan) — 전부 실제 임상 가이드라인(CPIC/FDA) 마커
- 질병 중심 재창출 후보 탐색: 간암 → PI3K/mTOR 억제제, CFH 연결 질병 →
  Pegcetacoplan/Iptacopan/Danicopan(실제 승인된 보체억제제) — 방법론이 실제 약물과 일치

## 아티팩트

- 📜 [Pathway to Pedigree](https://claude.ai/code/artifact/fbb4aa8c-875b-4c25-a8ac-3972394a80e1) — 전체 6단계 분석 보고서
- 🧬 [KEGG Relations Atlas](https://claude.ai/code/artifact/5dab17ec-13a4-4f0a-9203-ca4466ab5426) — 관계 그래프 허브 랭킹 + 인터랙티브 관계망
- 💊 [Drug Taxonomy](https://claude.ai/code/artifact/9305b69f-25a7-4c4f-a1b9-fa129b9bbb85) — ATC 통계 + dgroup 분류 트리

## 재현

```bash
cd "KEGG 데이터/script"
python build_kegg_tables.py
python analyze_kegg_relations.py
python disease_centered_analysis.py
python atc_classification_stats.py
python dgroup_class_tree.py

cd "../../가족 유전체(SNP) 데이터/script"
python family_snp_cross_analysis.py
# NCBI ClinVar variant_summary.txt.gz 다운로드 후 family_clinvar_match.py의 CLINVAR_GZ 경로 수정
python family_clinvar_match.py
python family_clinvar_kegg_crossvalidate.py
python family_pharmacogenomics.py
python family_kegg_gene_screening.py
python cfh_family_disease_profile.py
```

Python 3.11+, `pandas` 필요.

## 주의사항

이 프로젝트의 모든 "재창출 후보"·"동반질환"·"PGx 매칭"·"CFH 프로파일"은 **가설 생성용
스크리닝 결과**이며 의학적 진단·처방 근거가 아닙니다. 실제 임상 적용 전에는 반드시
의사·약사·유전상담사와 상의해야 합니다.
