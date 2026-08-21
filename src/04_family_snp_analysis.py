"""가족 SNP 데이터와 KEGG 질병-변이 데이터를 '안전한 기술 교차분석' 수준에서 연결한다.

이 스크립트는 두 가지만 한다:
  1. 가족 구성원별 SNP 원본 데이터의 기초 통계(콜 수, 결측률, 이형접합률, 염색체 분포)
  2. KEGG variant.txt에서 뽑은 909개 질병연관 유전자를 하나의 "관심 유전자 패널"로 정리

의도적으로 하지 않는 것: 특정 rsid가 특정 질병/변이와 일치하는지 판정하지 않는다.
KEGG variant.txt는 rsid가 아니라 OMIM 변이 번호 기준이고, 로컬 데이터에는 유전자의
genomic 좌표(어느 염색체 몇 번 위치인지)가 없어 SNP 좌표와 유전자를 연결할 수 없다.
즉 "이 사람이 이 유전자 영역에 SNP가 있다"조차 판정 불가 — 유전자 심볼 목록만 참고용으로 내보낸다.

출력:
  results/tables/family_snp_summary.csv       - 구성원별 SNP 기초 통계
  results/tables/family_snp_by_chromosome.csv - 구성원 x 염색체 SNP 개수
  data/processed/kegg_disease_gene_panel.csv  - KEGG variant.txt 기반 질병연관 유전자 패널(909개),
                                                 05_clinvar_annotation.py 등 뒤 단계에서 재사용됨
"""

import sys
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
FAMILY_DATA_DIR = ROOT / "data" / "raw" / "family_snp"
RESULTS_DIR = ROOT / "results" / "tables"
PROCESSED_DIR = ROOT / "data" / "processed"
KEGG_PROCESSED_DIR = PROCESSED_DIR

# data 폴더에는 "- Copy" 중복 파일, Child 3와 완전히 동일한 "Family Genome.csv",
# 가족 구성원이 아닌 공개 참조 게놈 "genome_zeeshan_usmani.csv"가 섞여 있어
# (diff로 바이트 단위까지 확인) 실제 가족 구성원 5명 파일만 명시적으로 지정한다.
FAMILY_MEMBER_FILES = {
    "Father": "Father Genome.csv",
    "Mother": "Mother Genome.csv",
    "Child 1": "Child 1 Genome.csv",
    "Child 2": "Child 2 Genome.csv",
    "Child 3": "Child 3 Genome.csv",
}

GENOME_COLUMNS = ["rsid", "chromosome", "position", "genotype"]


def load_genome(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig", comment="#", names=GENOME_COLUMNS, dtype=str)


def summarize_member(member: str, df: pd.DataFrame):
    genotype = df["genotype"].fillna("")
    length = genotype.str.len()

    no_call = genotype.isin(["--", "-"])
    diploid = (length == 2) & ~no_call
    haploid = (length == 1) & (genotype != "-")

    homozygous = diploid & (genotype.str[0] == genotype.str[1])
    heterozygous = diploid & (genotype.str[0] != genotype.str[1])

    diploid_calls = int(homozygous.sum() + heterozygous.sum())
    heterozygosity_rate = (heterozygous.sum() / diploid_calls) if diploid_calls else 0.0

    return {
        "member": member,
        "total_snps": len(df),
        "distinct_chromosomes": df["chromosome"].nunique(),
        "no_call_count": int(no_call.sum()),
        "no_call_rate": round(no_call.sum() / len(df), 5),
        "haploid_call_count": int(haploid.sum()),
        "homozygous_count": int(homozygous.sum()),
        "heterozygous_count": int(heterozygous.sum()),
        "heterozygosity_rate": round(heterozygosity_rate, 5),
    }


def build_family_snp_tables():
    summary_rows = []
    chrom_rows = []

    for member, filename in FAMILY_MEMBER_FILES.items():
        path = FAMILY_DATA_DIR / filename
        df = load_genome(path)

        summary_rows.append(summarize_member(member, df))

        chrom_counts = df["chromosome"].value_counts()
        for chrom, count in chrom_counts.items():
            chrom_rows.append({"member": member, "chromosome": chrom, "snp_count": int(count)})

        print(f"  {member:8s}: {len(df):,}행 ({filename})")

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(RESULTS_DIR / "family_snp_summary.csv", index=False, encoding="utf-8-sig")

    chrom_df = pd.DataFrame(chrom_rows)
    chrom_df.to_csv(RESULTS_DIR / "family_snp_by_chromosome.csv", index=False, encoding="utf-8-sig")

    return summary_df


def build_gene_panel():
    variant_df = pd.read_csv(KEGG_PROCESSED_DIR / "kegg_variant.csv", dtype=str).fillna("")
    disease_df = pd.read_csv(KEGG_PROCESSED_DIR / "kegg_disease.csv", dtype=str).fillna("")
    drug_df = pd.read_csv(KEGG_PROCESSED_DIR / "kegg_drug.csv", dtype=str).fillna("")

    disease_name = dict(zip(disease_df["entry_id"], disease_df["name"]))
    drug_name = dict(zip(drug_df["entry_id"], drug_df["name"]))

    def split_ids(cell: str):
        return [x for x in cell.split(";") if x]

    rows = []
    for gene_symbol, group in variant_df.groupby("gene_symbol"):
        if not gene_symbol:
            continue

        ko_ids = sorted({ko for cell in group["ko_id"] for ko in split_ids(cell)})

        disease_ids = sorted({d for cell in group["disease_ids"] for d in split_ids(cell)})
        drug_ids = sorted({d for cell in group["drug_target_ids"] for d in split_ids(cell)})
        network_ids = sorted({n for cell in group["network_ids"] for n in split_ids(cell)})

        disease_names = [disease_name.get(d, d).split(";")[0].strip() for d in disease_ids]
        drug_names = [drug_name.get(d, d).split(";")[0].strip() for d in drug_ids]

        rows.append(
            {
                "gene_symbol": gene_symbol,
                "ko_id": ";".join(ko_ids),
                "variant_entry_count": len(group),
                "variant_phenotypes": "; ".join(sorted(set(group["name"]))),
                "linked_disease_count": len(disease_ids),
                "linked_disease_names": "; ".join(disease_names),
                "linked_drug_target_count": len(drug_ids),
                "linked_drug_target_names": "; ".join(drug_names),
                "linked_network_count": len(network_ids),
            }
        )

    panel_df = pd.DataFrame(rows).sort_values(
        ["linked_disease_count", "variant_entry_count"], ascending=False
    )
    panel_df.to_csv(PROCESSED_DIR / "kegg_disease_gene_panel.csv", index=False, encoding="utf-8-sig")
    return panel_df


def main():
    print("가족 SNP 원본 통계:")
    summary_df = build_family_snp_tables()
    print()
    print(summary_df.to_string(index=False))
    print()

    panel_df = build_gene_panel()
    print(f"KEGG 질병연관 유전자 패널: {len(panel_df)}개 유전자")
    print(panel_df.head(10)[["gene_symbol", "variant_entry_count", "linked_disease_count", "linked_drug_target_count"]].to_string(index=False))
    print()
    print(f"-> {RESULTS_DIR / 'family_snp_summary.csv'}")
    print(f"-> {RESULTS_DIR / 'family_snp_by_chromosome.csv'}")
    print(f"-> {PROCESSED_DIR / 'kegg_disease_gene_panel.csv'}")
    print()
    print("주의: 이 패널은 '어떤 유전자가 KEGG에 질병연관 변이로 등록돼 있는지' 목록일 뿐,")
    print("가족 SNP 파일의 특정 rsid가 이 유전자들과 실제로 일치하는지는 판정하지 않았습니다.")
    print("(로컬 데이터에 유전자 genomic 좌표가 없어 rsid <-> 유전자 매칭이 불가능합니다.)")


if __name__ == "__main__":
    main()
