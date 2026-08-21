"""data/raw/kegg 안의 5개 flat file(dgroup, disease, drug, network, variant)을
구조화된 CSV 테이블로 변환한다.

출력 (data/processed/):
  kegg_dgroup.csv    - 약물군(DGroup)
  kegg_disease.csv   - 질병(Disease)
  kegg_drug.csv      - 개별 약물(Drug)
  kegg_network.csv   - 신호전달/질병 네트워크(Network)
  kegg_variant.csv   - 유전자 변이(Variant)
  kegg_relations.csv - 엔트리 간 관계 edge list (source -> target)
                        network/그래프 분석(NetworkX, Gephi 등)에 바로 사용 가능
"""

import csv
import re
import sys
from pathlib import Path

if sys.stdout.encoding and sys.stdout.encoding.lower() != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

from kegg_flatfile_parser import parse_entries, entry_id_and_type, field_text, field_lines

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "raw" / "kegg"
OUTPUT_DIR = ROOT / "data" / "processed"

RE_H = re.compile(r"\bH\d{5}\b")
RE_D = re.compile(r"\bD\d{5}\b")
RE_DG = re.compile(r"\bDG\d{5}\b")
RE_N = re.compile(r"\bN\d{5}\b")
RE_C = re.compile(r"\bC\d{5}\b")
RE_K = re.compile(r"\bK\d{5}\b")
RE_VARIANT = re.compile(r"\b\d+v\d+\b")
RE_BRACKET = re.compile(r"\[([A-Z]+):([^\]]+)\]")


def uniq(seq):
    seen = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def ids_by_pattern(text: str, pattern: re.Pattern):
    return uniq(pattern.findall(text))


def id_type(id_str: str) -> str:
    if RE_DG.fullmatch(id_str):
        return "dgroup"
    if RE_D.fullmatch(id_str):
        return "drug"
    if RE_H.fullmatch(id_str):
        return "disease"
    if RE_N.fullmatch(id_str):
        return "network"
    if RE_C.fullmatch(id_str):
        return "compound"
    if RE_K.fullmatch(id_str):
        return "ortholog"
    if RE_VARIANT.fullmatch(id_str):
        return "variant"
    if id_str.isdigit():
        return "gene_entrez"
    return "unknown"


def write_csv(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


class RelationCollector:
    def __init__(self):
        self.edges = []

    def add(self, src_type, src_id, dst_id, field):
        dst_type = id_type(dst_id)
        if dst_type == "unknown" or dst_id == src_id:
            return
        self.edges.append((src_type, src_id, dst_type, dst_id, field))

    def add_many(self, src_type, src_id, dst_ids, field):
        for dst_id in dst_ids:
            self.add(src_type, src_id, dst_id, field)

    def write(self, path: Path):
        rows = uniq(self.edges)
        write_csv(
            path,
            ["source_type", "source_id", "target_type", "target_id", "field"],
            rows,
        )
        return len(rows)


def parse_key_value_lines(lines) -> dict:
    """'Key: value' 형태로 된 줄들(REMARK, DBLINKS 등)을 dict로 변환한다.
    각 줄이 하나의 key: value 쌍이라는 KEGG flat file 관례에 따라 줄 단위로 분리한다."""
    result = {}
    for line in lines:
        if ":" in line:
            key, _, value = line.partition(":")
            key = key.strip()
            value = value.strip()
            if key and key not in result:
                result[key] = value
    return result


def build_dgroup(rel: RelationCollector):
    rows = []
    for entry in parse_entries(DATA_DIR / "dgroup.txt"):
        entry_id, entry_type = entry_id_and_type(entry)
        if not entry_id:
            continue
        name = field_text(entry, "NAME")
        dtype = field_text(entry, "TYPE")
        class_text = field_text(entry, "CLASS")
        remark = parse_key_value_lines(field_lines(entry, "REMARK"))
        comment = field_text(entry, "COMMENT")

        member_ids = ids_by_pattern(field_text(entry, "MEMBER"), RE_D)
        class_dg_ids = ids_by_pattern(class_text, RE_DG)
        atc_codes = remark.get("ATC code", "")

        rel.add_many("dgroup", entry_id, member_ids, "MEMBER")
        rel.add_many("dgroup", entry_id, class_dg_ids, "CLASS")

        rows.append(
            [entry_id, name, dtype, class_text, len(member_ids), ";".join(member_ids), atc_codes, comment]
        )

    write_csv(
        OUTPUT_DIR / "kegg_dgroup.csv",
        ["entry_id", "name", "type", "class", "member_count", "member_ids", "atc_codes", "comment"],
        rows,
    )
    return len(rows)


def build_disease(rel: RelationCollector):
    rows = []
    for entry in parse_entries(DATA_DIR / "disease.txt"):
        entry_id, entry_type = entry_id_and_type(entry)
        if not entry_id:
            continue
        name = field_text(entry, "NAME", sep=" ")
        category = field_text(entry, "CATEGORY")
        description = field_text(entry, "DESCRIPTION")
        dblinks = field_text(entry, "DBLINKS")
        pathogen = field_text(entry, "PATHOGEN")

        gene_lines = field_lines(entry, "GENE")
        gene_symbols = uniq(line.split()[0] for line in gene_lines if line.split())
        gene_hsa_ids = ids_by_pattern(field_text(entry, "GENE"), re.compile(r"\[HSA:(\d+)\]"))

        drug_ids = ids_by_pattern(field_text(entry, "DRUG"), re.compile(r"\[DR:(D\d{5})\]"))
        network_ids = ids_by_pattern(field_text(entry, "NETWORK"), RE_N)

        rel.add_many("disease", entry_id, gene_hsa_ids, "GENE")
        rel.add_many("disease", entry_id, drug_ids, "DRUG")
        rel.add_many("disease", entry_id, network_ids, "NETWORK")

        rows.append(
            [
                entry_id,
                name,
                category,
                description,
                len(gene_symbols),
                ";".join(gene_symbols),
                len(drug_ids),
                ";".join(drug_ids),
                ";".join(network_ids),
                dblinks,
                pathogen,
            ]
        )

    write_csv(
        OUTPUT_DIR / "kegg_disease.csv",
        [
            "entry_id",
            "name",
            "category",
            "description",
            "gene_count",
            "gene_symbols",
            "drug_count",
            "drug_ids",
            "network_ids",
            "dblinks",
            "pathogen",
        ],
        rows,
    )
    return len(rows)


def build_drug(rel: RelationCollector):
    rows = []
    for entry in parse_entries(DATA_DIR / "drug .txt"):
        entry_id, entry_type = entry_id_and_type(entry)
        if not entry_id:
            continue
        name = field_text(entry, "NAME", sep=" ")
        formula = field_text(entry, "FORMULA")
        mol_weight = field_text(entry, "MOL_WEIGHT")
        efficacy = field_text(entry, "EFFICACY")
        target_text = field_text(entry, "TARGET")
        remark = parse_key_value_lines(field_lines(entry, "REMARK"))
        dblinks = parse_key_value_lines(field_lines(entry, "DBLINKS"))

        atc_codes = remark.get("ATC code", "")
        dgroup_id = remark.get("Chemical structure group", "")
        compound_id = remark.get("Same as", "")
        therapeutic_category = remark.get("Therapeutic category", "")

        target_hsa_ids = ids_by_pattern(target_text, re.compile(r"\[HSA:([\d\s]+)\]"))
        target_hsa_ids = uniq(
            gid for group in target_hsa_ids for gid in group.split()
        )

        if dgroup_id:
            rel.add("drug", entry_id, dgroup_id, "REMARK:Chemical structure group")
        if compound_id:
            rel.add("drug", entry_id, compound_id, "REMARK:Same as")
        rel.add_many("drug", entry_id, target_hsa_ids, "TARGET")

        rows.append(
            [
                entry_id,
                name,
                formula,
                mol_weight,
                atc_codes,
                dgroup_id,
                compound_id,
                therapeutic_category,
                efficacy,
                dblinks.get("CAS", ""),
                dblinks.get("PubChem", ""),
                dblinks.get("ChEBI", ""),
            ]
        )

    write_csv(
        OUTPUT_DIR / "kegg_drug.csv",
        [
            "entry_id",
            "name",
            "formula",
            "mol_weight",
            "atc_codes",
            "dgroup_id",
            "compound_id",
            "therapeutic_category",
            "efficacy",
            "cas",
            "pubchem",
            "chebi",
        ],
        rows,
    )
    return len(rows)


def build_network(rel: RelationCollector):
    rows = []
    for entry in parse_entries(DATA_DIR / "network.txt"):
        entry_id, entry_type = entry_id_and_type(entry)
        if not entry_id:
            continue
        name = field_text(entry, "NAME")
        ntype = field_text(entry, "TYPE")
        definition = field_text(entry, "DEFINITION")
        class_lines = field_lines(entry, "CLASS")

        gene_lines = field_lines(entry, "GENE")
        gene_entrez_ids = uniq(line.split()[0] for line in gene_lines if line.split())
        gene_symbols = uniq(
            line.split(None, 1)[1].split(";")[0].strip()
            for line in gene_lines
            if len(line.split(None, 1)) > 1
        )

        disease_ids = ids_by_pattern(field_text(entry, "DISEASE"), RE_H)
        variant_ids = uniq(line.split()[0] for line in field_lines(entry, "VARIANT") if line.split())

        rel.add_many("network", entry_id, gene_entrez_ids, "GENE")
        rel.add_many("network", entry_id, disease_ids, "DISEASE")
        rel.add_many("network", entry_id, variant_ids, "VARIANT")

        rows.append(
            [
                entry_id,
                name,
                ntype,
                definition,
                "; ".join(class_lines),
                len(gene_entrez_ids),
                ";".join(gene_symbols),
                ";".join(disease_ids),
                ";".join(variant_ids),
            ]
        )

    write_csv(
        OUTPUT_DIR / "kegg_network.csv",
        [
            "entry_id",
            "name",
            "type",
            "definition",
            "class",
            "gene_count",
            "gene_symbols",
            "disease_ids",
            "variant_ids",
        ],
        rows,
    )
    return len(rows)


def build_variant(rel: RelationCollector):
    rows = []
    for entry in parse_entries(DATA_DIR / "variant.txt"):
        entry_id, entry_type = entry_id_and_type(entry)
        if not entry_id:
            continue
        name = field_text(entry, "NAME")
        organism = field_text(entry, "ORGANISM")
        variation = field_text(entry, "VARIATION")

        gene_line = field_text(entry, "GENE")
        gene_symbol = gene_line.split()[0] if gene_line else ""
        ko_ids = ids_by_pattern(gene_line, RE_K)

        disease_ids = ids_by_pattern(field_text(entry, "DISEASE"), RE_H)
        network_ids = ids_by_pattern(field_text(entry, "NETWORK"), RE_N)
        drug_target_ids = ids_by_pattern(field_text(entry, "DRUG_TARGET"), RE_D)

        rel.add_many("variant", entry_id, ko_ids, "GENE")
        rel.add_many("variant", entry_id, disease_ids, "DISEASE")
        rel.add_many("variant", entry_id, network_ids, "NETWORK")
        rel.add_many("variant", entry_id, drug_target_ids, "DRUG_TARGET")

        rows.append(
            [
                entry_id,
                name,
                gene_symbol,
                ";".join(ko_ids),
                organism,
                variation,
                ";".join(disease_ids),
                ";".join(network_ids),
                ";".join(drug_target_ids),
            ]
        )

    write_csv(
        OUTPUT_DIR / "kegg_variant.csv",
        [
            "entry_id",
            "name",
            "gene_symbol",
            "ko_id",
            "organism",
            "variation",
            "disease_ids",
            "network_ids",
            "drug_target_ids",
        ],
        rows,
    )
    return len(rows)


def main():
    rel = RelationCollector()

    counts = {
        "dgroup": build_dgroup(rel),
        "disease": build_disease(rel),
        "drug": build_drug(rel),
        "network": build_network(rel),
        "variant": build_variant(rel),
    }

    edge_count = rel.write(OUTPUT_DIR / "kegg_relations.csv")

    print("엔트리 파싱 결과:")
    for name, count in counts.items():
        print(f"  {name}: {count}건")
    print(f"관계(edge): {edge_count}건")
    print(f"-> {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
