"""variant.txt(KEGG variant flat file)에서 유전자 목록을 추출한다.

입력: ../data/variant.txt
출력:
  ../output/variant_genes.csv         - 변이 엔트리별 유전자 매핑 (1행 = 1 variant entry)
  ../output/variant_genes_unique.csv  - 중복 제거된 유전자 목록 (유전자별 관련 variant entry 수)
"""

import csv
import re
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR.parent / "data" / "variant.txt"
OUTPUT_DIR = SCRIPT_DIR.parent / "output"

ENTRY_RE = re.compile(r"^ENTRY\s+(\S+)")
GENE_RE = re.compile(r"^GENE\s+(\S+)\s+(.*)$")
KO_RE = re.compile(r"\[KO:([^\]]+)\]")


def parse_variant_file(path: Path):
    """variant.txt를 한 줄씩 읽어 (variant_id, gene_symbol, gene_description, ko_ids) 리스트 반환."""
    records = []
    entry_id = None

    with open(path, encoding="utf-8") as f:
        for line in f:
            entry_match = ENTRY_RE.match(line)
            if entry_match:
                entry_id = entry_match.group(1)
                continue

            gene_match = GENE_RE.match(line)
            if gene_match and entry_id:
                symbol = gene_match.group(1)
                rest = gene_match.group(2)

                ko_match = KO_RE.search(rest)
                ko_ids = ko_match.group(1) if ko_match else ""
                description = KO_RE.sub("", rest).strip()

                records.append((entry_id, symbol, description, ko_ids))

    return records


def write_per_entry_csv(records, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["variant_id", "gene_symbol", "gene_description", "ko_id"])
        writer.writerows(records)


def write_unique_gene_csv(records, path: Path):
    counts = {}
    for _, symbol, description, ko_ids in records:
        if symbol not in counts:
            counts[symbol] = {"description": description, "ko_id": ko_ids, "count": 0}
        counts[symbol]["count"] += 1

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_symbol", "gene_description", "ko_id", "variant_entry_count"])
        for symbol in sorted(counts):
            info = counts[symbol]
            writer.writerow([symbol, info["description"], info["ko_id"], info["count"]])

    return len(counts)


def main():
    records = parse_variant_file(INPUT_FILE)

    per_entry_path = OUTPUT_DIR / "variant_genes.csv"
    unique_path = OUTPUT_DIR / "variant_genes_unique.csv"

    write_per_entry_csv(records, per_entry_path)
    unique_count = write_unique_gene_csv(records, unique_path)

    print(f"variant entries parsed: {len(records)}")
    print(f"unique genes: {unique_count}")
    print(f"-> {per_entry_path}")
    print(f"-> {unique_path}")


if __name__ == "__main__":
    main()
