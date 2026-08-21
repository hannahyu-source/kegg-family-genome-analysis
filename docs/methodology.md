# Methodology

이 문서는 `src/`의 각 파이프라인 스크립트가 실제로 무엇을 하는지 설명한다. 실행 순서와
입출력 파일 목록은 [`workflow.md`](workflow.md)를 참고하고, 각 단계에서 발견한 버그·구체적
수치·판단 근거는 [`development_log.md`](development_log.md)를 참고할 것.

## 1. KEGG 데이터 파싱

**스크립트**: `src/kegg_flatfile_parser.py`(공용 파서), `src/01_build_kegg_tables.py`

KEGG flat file은 필드명이 앞 12칸에 오고 내용이 13번째 칸부터 시작하는 고정폭 텍스트
포맷이며, 엔트리는 `///` 줄로 구분된다. 필드명이 비어 있는 줄(앞 12칸이 공백)은 바로 위
필드의 연속 줄로 취급한다.

`dgroup.txt` / `disease.txt` / `drug .txt` / `network.txt` / `variant.txt` 5개 원본 파일에서
각각 다음 필드를 추출한다:

- **dgroup**: NAME, TYPE, CLASS(상위 분류), MEMBER(소속 drug ID), REMARK의 ATC code
- **disease**: NAME, CATEGORY, GENE([HSA:entrez_id] 형태), DRUG([DR:D#####]), NETWORK(N#####)
- **drug**: NAME, FORMULA, ATC code, TARGET([HSA:...]), REMARK의 dgroup/compound 참조,
  DBLINKS(CAS/PubChem/ChEBI)
- **network**: NAME, TYPE, GENE(Entrez ID + 심볼), DISEASE(H#####), VARIANT
- **variant**: NAME, GENE 심볼, KO ID, DISEASE(H#####), NETWORK(N#####), DRUG_TARGET(D#####)

각 필드 안에서 정규식으로 다른 엔트리 타입의 ID(`H#####`=disease, `D#####`=drug,
`DG#####`=dgroup, `N#####`=network, `C#####`=compound, `K#####`=ortholog, 숫자만=Entrez 유전자
ID)를 모두 추출해 **source_type, source_id, target_type, target_id, field** 5개 컬럼의 edge
list(`kegg_relations.csv`)로 통합한다. 이 관계는 KEGG가 명시한 방향(예: disease→gene,
network→disease)을 그대로 따르며, 반대 방향으로 추론하지 않는다.

**엔트리 타입별 개수**(2026-08 다운로드 시점): dgroup 2,429 / disease 2,633 / drug 12,115 /
network 1,417 / variant 982, 관계 69,546건.

## 2. 가족 SNP 전처리

**스크립트**: `src/04_family_snp_analysis.py`

23andMe 원시 genotype CSV는 `rsid, chromosome, position, genotype` 4컬럼, 헤더 없이 `#`
주석으로 시작한다. 원본 데이터셋의 13개 파일 중 실제 가족 구성원은 5명(Father/Mother/
Child 1/2/3)뿐이며, 나머지(`*-Copy.csv` 중복, Child 3와 바이트 단위로 동일한
`Family Genome.csv`, 가족이 아닌 공개 참조 게놈 `genome_zeeshan_usmani.csv`)는 diff로
확인 후 명시적으로 제외했다.

기초 통계는 다음과 같이 정의한다:
- **무호출(no-call)**: genotype이 `--` 또는 `-`
- **이형접합(heterozygous)**: 2글자 genotype이고 두 글자가 다름
- **이형접합률**: 이형접합 개수 / (동형접합 + 이형접합 개수), 무호출 제외

결측값은 별도 대체(impute)하지 않고 `no_call_count`/`no_call_rate`로 집계만 한다. 가족
구성원 간 병합은 이 단계에서 하지 않으며(각자 독립적으로 통계만 낸다), 실제 병합은 3단계
ClinVar 매칭에서 rsid를 공통 키로 삼아 이루어진다.

## 3. ClinVar 매칭

**스크립트**: `src/05_clinvar_annotation.py`

- **소스 파일**: NCBI ClinVar `variant_summary.txt` (탭 구분, 약 900만 행)
- **genome build**: GRCh37만 사용(`Assembly == "GRCh37"`) — 가족 23andMe 데이터의 기준
  빌드와 맞추기 위함
- **매칭 키**: `RS# (dbSNP)` 컬럼을 `rs` 접두사를 붙여 `rsid`로 변환한 뒤, 가족 5명의 rsid
  합집합(1,117,586개)에 속하는 행만 유지
- **필터**: RS#가 없는 행(`-1`)은 제외. 대용량 파일이라 20만 행 단위 청크로 스트리밍 필터링
- **allele 비교**: SNV(단일 염기 변이)이고 Ref/Alt가 각각 한 글자인 경우에 한해서만
  "이 사람의 genotype 문자열에 ClinVar의 alternate allele 문자가 포함돼 있는가
  (`carries_alt_allele`)"를 계산한다. indel·복합 변이는 23andMe의 D/I 표기와 VCF 스타일
  Ref/Alt를 직접 비교할 근거가 없어 계산하지 않고 `None`으로 남긴다.

**rsid 수준 매칭의 한계**: `carries_alt_allele=True`는 "그 사람의 2개 대립유전자 문자열 안에
ClinVar가 등록한 대체 대립유전자 문자가 존재한다"는 뜻이지, 염색체 phase(부모 중 어느 쪽에서
왔는지)나 정확한 zygosity(동형/이형 여부)를 항상 명확히 구분하지 않는다. 또한 이 매칭은
rsid 단위이며 chromosome/position/REF/ALT를 이용한 완전한 좌표 정규화는 하지 않는다 — 자세한
한계는 [`limitations.md`](limitations.md)를 참고.

## 4. KEGG 교차검증

**스크립트**: `src/06_kegg_crossvalidation.py`

ClinVar에서 병원성 계열(`Pathogenic` 포함, `Conflicting` 제외)이고 가족 중 누군가 실제로
대체 대립유전자를 보유한 것으로 확인된 변이만 후보로 남긴다. 이 유전자 심볼이
`kegg_disease_gene_panel.csv`(KEGG `variant.txt`가 실제로 `DISEASE: H#####` 필드를 명시한
909개 유전자 패널)에도 있으면 "두 개의 독립적으로 큐레이션된 데이터베이스가 모두 이
유전자를 질병과 연결짓고 있다"는 의미로 교차검증 통과 처리한다.

이 교차검증이 확인하는 것은 **"같은 유전자가 두 DB에서 질병과 연결된다"**는 것이지,
**"ClinVar가 병원성으로 판정한 바로 그 변이의 기전을 KEGG가 확인해준다"**는 뜻이 아니다 —
유전자 수준의 독립 근거이지, 변이 수준의 인과 증명이 아니다.

## 5. 가족 비교

**스크립트**: `src/04_family_snp_analysis.py`(기초 통계), `src/08_gene_screening.py`(유전자
단위 보유 여부)

가족 비교는 다음 두 수준에서만 이루어진다:
- **구성원별 독립 통계**: 총 SNP 수, 결측률, 이형접합률 (4단계)
- **rsid/유전자 단위 보유자 비교**: 각 ClinVar 매칭 변이에 대해 `{member}_carries_alt_allele`
  플래그로 어느 구성원이 그 대립유전자를 보유했는지 나열(`carriers` 컬럼)

**하지 않는 것**: 이 프로젝트는 부모-자식 유전 패턴에 대한 정식 멘델 분리 분석(Mendelian
segregation analysis)이나 친자 확인, 좌표 기반 haplotype phasing을 수행하지 않는다. "Father와
Child 1이 같은 rsid를 보유"는 유전됐을 가능성을 시사할 뿐, 그 자체로 유전 경로를 증명하지
않는다.

## 6. 약물유전체(PGx)

**스크립트**: `src/07_pharmacogenomics.py`

ClinVar `ClinicalSignificance`가 `"drug response"`를 포함하는 변이 중 가족이 실제로 보유한
것만 추출한다(60건, 유전자 26개). KEGG와는 두 가지 방식으로 연결한다:

1. **유전자 심볼 → KEGG dgroup 이름 매칭**: 예를 들어 유전자가 `CYP2D6`이면 이름에
   "CYP2D6"가 들어간 모든 dgroup(예: "CYP2D6 substrate")과 그 dgroup에 속한 모든 약물을
   넓은 후보군으로 잡는다.
2. **PhenotypeList → KEGG drug 이름 직접 매칭**: ClinVar의 `PhenotypeList`에 있는
   `"약물명 response"` 패턴(정규식 `[A-Za-z][A-Za-z0-9/\-]{2,}\s+response`)에서 약물명을
   뽑아 KEGG `kegg_drug.csv`의 대표 이름(괄호/세미콜론 이전 부분)과 직접 매칭한다. 이 방식이
   훨씬 좁고 구체적인 연결이다(53/60건 성공).

## 7. 사례 연구(case study) 선정 근거

**스크립트**: `src/09_cfh_case_study.py`

CFH를 사례 연구로 선정한 이유는 사전에 정한 것이 아니라, 3~6단계 파이프라인이 실제로
찾아낸 결과이기 때문이다: 가족이 실제로 보유한 ClinVar 병원성(비상충) 변이 4건(CFH×2,
LMNA×1, UGT1A1×1) 중, **KEGG `variant.txt`에도 명시적 질병 링크가 있어 교차검증을 통과한
유전자는 CFH뿐**이었다(LMNA·UGT1A1은 KEGG 쪽에 해당 병원성 맥락의 명시적 질병 링크가 없어
탈락). 즉 CFH는 "흥미로워 보여서" 고른 것이 아니라, 이 프로젝트가 정의한 교차검증 기준을
통과한 유일한 유전자라서 선정했다.
