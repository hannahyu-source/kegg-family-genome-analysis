"""kegg_drug.csv의 ATC 코드로 약물을 치료 영역별로 집계한다.

ATC(Anatomical Therapeutic Chemical) 코드는 WHO 표준 분류로, 첫 글자가
해부학적 주분류(예: A=소화기, C=심혈관계, N=신경계)를 나타낸다. 한 약물이
여러 ATC 코드(=여러 치료 영역)를 가질 수 있어, "약물 1개당 중복 집계"와
"영역별 고유 약물 수" 두 가지를 모두 낸다.

출력 (results/tables/):
  atc_top_level_summary.csv  - 1단계(해부학적 주분류) 약물 수
  atc_subgroup_summary.csv   - 3단계(치료학적 소분류) 약물 수 상위 30
"""

import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUT_DIR = ROOT / "results" / "tables"

ATC_TOP_LEVEL = {
    "A": "소화기·대사(Alimentary tract & metabolism)",
    "B": "혈액·조혈기관 (Blood & blood forming organs)",
    "C": "심혈관계 (Cardiovascular system)",
    "D": "피부 (Dermatologicals)",
    "G": "비뇨생식기·성호르몬 (Genito-urinary system & sex hormones)",
    "H": "전신 호르몬제 (Systemic hormonal preparations)",
    "J": "전신 항감염제 (Antiinfectives for systemic use)",
    "L": "항종양·면역조절제 (Antineoplastic & immunomodulating agents)",
    "M": "근골격계 (Musculo-skeletal system)",
    "N": "신경계 (Nervous system)",
    "P": "항기생충제 (Antiparasitic products)",
    "R": "호흡기계 (Respiratory system)",
    "S": "감각기관 (Sensory organs)",
    "V": "기타 (Various)",
}

ATC_CODE_RE = re.compile(r"\b([A-Z])(\d{2})([A-Z]{2})(\d{2})\b")


def main():
    drug_df = pd.read_csv(PROCESSED_DIR / "kegg_drug.csv", dtype=str).fillna("")
    with_atc = drug_df[drug_df["atc_codes"] != ""]
    print(f"ATC 코드가 있는 약물: {len(with_atc)} / {len(drug_df)}")

    top_drug_sets = {k: set() for k in ATC_TOP_LEVEL}
    top_code_counter = Counter()
    sub_drug_sets = {}

    for entry_id, codes in zip(with_atc["entry_id"], with_atc["atc_codes"]):
        for code in codes.split():
            m = ATC_CODE_RE.match(code)
            if not m:
                continue
            top = m.group(1)
            sub3 = f"{m.group(1)}{m.group(2)}{m.group(3)}"  # 3단계: 치료학적 소분류
            top_code_counter[top] += 1
            if top in top_drug_sets:
                top_drug_sets[top].add(entry_id)
            sub_drug_sets.setdefault(sub3, set()).add(entry_id)

    top_rows = [
        {
            "atc_top": k,
            "category": ATC_TOP_LEVEL[k],
            "unique_drug_count": len(top_drug_sets[k]),
            "atc_code_instances": top_code_counter[k],
        }
        for k in ATC_TOP_LEVEL
    ]
    top_df = pd.DataFrame(top_rows).sort_values("unique_drug_count", ascending=False)
    top_df.to_csv(OUTPUT_DIR / "atc_top_level_summary.csv", index=False, encoding="utf-8-sig")

    sub_rows = [
        {"atc_subgroup": k, "unique_drug_count": len(v)}
        for k, v in sub_drug_sets.items()
    ]
    sub_df = pd.DataFrame(sub_rows).sort_values("unique_drug_count", ascending=False).head(30)
    sub_df.to_csv(OUTPUT_DIR / "atc_subgroup_summary.csv", index=False, encoding="utf-8-sig")

    print()
    print(top_df.to_string(index=False))
    print()
    print("=== 3단계 소분류 상위 10 ===")
    print(sub_df.head(10).to_string(index=False))
    print()
    print(f"-> {OUTPUT_DIR / 'atc_top_level_summary.csv'}")
    print(f"-> {OUTPUT_DIR / 'atc_subgroup_summary.csv'}")


if __name__ == "__main__":
    main()
