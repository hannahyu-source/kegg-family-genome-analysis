"""가족 SNP 데이터를 NCBI ClinVar variant_summary.txt와 실제 rsid 단위로 매칭한다.

이전 단계(family_snp_cross_analysis.py)는 KEGG variant.txt에 rsid/좌표가 없어
"유전자 패널만 참고용으로 정리"하는 수준에서 멈췄다. ClinVar의 variant_summary.txt는
RS# (dbSNP), Assembly, GeneSymbol, ClinicalSignificance, ReferenceAlleleVCF/
AlternateAlleleVCF를 갖고 있어 실제 join이 가능하다.

방법:
  1. GRCh37(가족 23andMe 데이터의 기준 빌드) 행만 사용
  2. RS# (dbSNP)가 있고, 그 rsid가 가족 SNP 파일 어딘가에 존재하는 행만 유지
  3. SNV(단일 염기 변이)이고 Ref/Alt가 한 글자씩인 경우에 한해서만
     "이 사람 genotype에 ClinVar의 alternate allele이 포함돼 있는가"를 계산
     (indel/복합변이는 23andMe genotype 표기(D/I)와 직접 비교할 근거가 없어 계산하지 않음)

9백만 행짜리 파일이라 청크 단위로 스트리밍 필터링한다.

출력 (가족 유전체(SNP) 데이터/output/):
  family_clinvar_matches.csv - rsid 단위 매칭 결과, 구성원별 genotype + carries_alt_allele
  family_clinvar_summary.csv - 구성원별 · 임상적 유의성별 매칭 건수 요약
"""

import gzip
import pickle
import sys
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
FAMILY_OUTPUT_DIR = SCRIPT_DIR.parent / "output"
CLINVAR_GZ = Path(
    r"C:\Users\yuhan\AppData\Local\Temp\claude\C--A-track-hannah-KEGG-genomes-networks-diseases-and-drugs"
    r"\f1274abf-9e23-4eef-b0f4-82d6829cb3ef\scratchpad\clinvar\variant_summary.txt.gz"
)

MEMBER_ORDER = ["Father", "Mother", "Child 1", "Child 2", "Child 3"]

USECOLS = [
    "RS# (dbSNP)",
    "Assembly",
    "Type",
    "GeneSymbol",
    "ClinicalSignificance",
    "ClinSigSimple",
    "ReviewStatus",
    "NumberSubmitters",
    "PhenotypeList",
    "VariationID",
    "Chromosome",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
]

CHUNK_SIZE = 200_000


def load_family_lookup():
    with open(FAMILY_OUTPUT_DIR / "_family_genotypes.pkl", "rb") as f:
        per_member = pickle.load(f)
    with open(FAMILY_OUTPUT_DIR / "_family_rsid_set.pkl", "rb") as f:
        rsid_set = pickle.load(f)
    return per_member, rsid_set


def stream_matches(rsid_set: set):
    matched_chunks = []
    total_rows = 0
    total_grch37 = 0

    reader = pd.read_csv(
        CLINVAR_GZ,
        sep="\t",
        usecols=USECOLS,
        dtype=str,
        chunksize=CHUNK_SIZE,
        compression="gzip",
        na_filter=False,
    )

    for i, chunk in enumerate(reader):
        total_rows += len(chunk)
        chunk = chunk[chunk["Assembly"] == "GRCh37"]
        total_grch37 += len(chunk)

        rs = chunk["RS# (dbSNP)"]
        has_rs = rs != "-1"
        chunk = chunk[has_rs]
        if len(chunk) == 0:
            continue

        chunk = chunk.assign(rsid=("rs" + chunk["RS# (dbSNP)"]))
        chunk = chunk[chunk["rsid"].isin(rsid_set)]

        if len(chunk):
            matched_chunks.append(chunk)

        if (i + 1) % 5 == 0:
            print(f"  ...{total_rows:,}행 처리, 지금까지 매칭 {sum(len(c) for c in matched_chunks):,}건")

    print(f"전체 {total_rows:,}행 중 GRCh37 {total_grch37:,}행, 가족 rsid와 매칭 {sum(len(c) for c in matched_chunks):,}건")
    return pd.concat(matched_chunks, ignore_index=True) if matched_chunks else pd.DataFrame(columns=USECOLS + ["rsid"])


def add_family_genotypes(matches: pd.DataFrame, per_member: dict) -> pd.DataFrame:
    for member in MEMBER_ORDER:
        lookup = per_member[member]
        matches[f"{member}_genotype"] = matches["rsid"].map(lookup)
    return matches


def add_allele_check(matches: pd.DataFrame) -> pd.DataFrame:
    is_snv = matches["Type"].str.lower().str.contains("single nucleotide")
    ref = matches["ReferenceAlleleVCF"]
    alt = matches["AlternateAlleleVCF"]
    comparable = is_snv & (ref.str.len() == 1) & (alt.str.len() == 1)
    matches["allele_comparable"] = comparable

    for member in MEMBER_ORDER:
        genotype_col = f"{member}_genotype"
        flag_col = f"{member}_carries_alt_allele"
        matches[flag_col] = [
            (c and isinstance(g, str) and a in g) if c else None
            for c, g, a in zip(comparable, matches[genotype_col], alt)
        ]
    return matches


def main():
    print("가족 rsid 목록 / genotype 로딩 중...")
    per_member, rsid_set = load_family_lookup()
    print(f"  가족 전체 고유 rsid: {len(rsid_set):,}개")

    cache_path = FAMILY_OUTPUT_DIR / "_clinvar_raw_matches.pkl"
    if cache_path.exists():
        print("이전에 필터링해둔 ClinVar 매칭 결과 캐시를 재사용합니다...")
        matches = pd.read_pickle(cache_path)
    else:
        print("ClinVar variant_summary.txt 스트리밍 필터링 중 (9,044,811행)...")
        matches = stream_matches(rsid_set)
        matches.to_pickle(cache_path)

    print("가족 genotype 결합 중...")
    matches = add_family_genotypes(matches, per_member)

    print("allele 비교(SNV만) 계산 중...")
    matches = add_allele_check(matches)

    out_cols = [
        "rsid", "GeneSymbol", "Chromosome", "PositionVCF",
        "ReferenceAlleleVCF", "AlternateAlleleVCF", "Type",
        "ClinicalSignificance", "ClinSigSimple", "ReviewStatus", "NumberSubmitters",
        "PhenotypeList", "VariationID", "allele_comparable",
    ] + [f"{m}_genotype" for m in MEMBER_ORDER] + [f"{m}_carries_alt_allele" for m in MEMBER_ORDER]

    matches = matches[out_cols].sort_values(["ClinSigSimple", "NumberSubmitters"], ascending=[False, False])
    matches.to_csv(FAMILY_OUTPUT_DIR / "family_clinvar_matches.csv", index=False, encoding="utf-8-sig")

    # 요약: 구성원 x 임상적 유의성별 "alt allele 보유" 건수
    summary_rows = []
    for member in MEMBER_ORDER:
        flag_col = f"{member}_carries_alt_allele"
        carrying = matches[matches[flag_col] == True]
        for sig, group in carrying.groupby("ClinicalSignificance"):
            summary_rows.append({"member": member, "clinical_significance": sig, "variant_count": len(group)})
    summary_df = pd.DataFrame(summary_rows).sort_values(["member", "variant_count"], ascending=[True, False])
    summary_df.to_csv(FAMILY_OUTPUT_DIR / "family_clinvar_summary.csv", index=False, encoding="utf-8-sig")

    print()
    print(f"매칭 rsid 총 {len(matches):,}건 (allele 비교 가능: {int(matches['allele_comparable'].sum()):,}건)")
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_clinvar_matches.csv'}")
    print(f"-> {FAMILY_OUTPUT_DIR / 'family_clinvar_summary.csv'}")


if __name__ == "__main__":
    main()
