"""KEGG 분석 스크립트 여러 개가 공유하는 유틸리티 함수.

파이프라인 스크립트 파일명이 실행 순서를 나타내는 숫자로 시작해서
(예: 02_analyze_kegg_relations.py) 서로를 `import` 문으로 직접 불러올 수
없다 (Python 문법상 식별자가 숫자로 시작할 수 없음). 그래서 여러 스크립트가
공유하는 함수는 숫자 접두사가 없는 이 모듈에 모아두고, 각 파이프라인
스크립트가 여기서 import한다.
"""

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from kegg_flatfile_parser import parse_entries


def build_gene_symbol_map(data_dir: Path) -> dict:
    """network.txt/disease.txt의 GENE 필드에서 Entrez ID -> 유전자 심볼 매핑을 뽑는다.
    두 파일 모두 'entrez_id  SYMBOL; description' 또는 'SYMBOL ... [HSA:entrez_id]' 형태로 심볼을 담고 있다."""
    mapping = {}

    for entry in parse_entries(data_dir / "network.txt"):
        for line in entry.get("GENE", []):
            parts = line.strip().split(None, 1)
            if len(parts) == 2 and parts[0].isdigit():
                mapping.setdefault(parts[0], parts[1].split(";")[0].strip())

    hsa_pattern = re.compile(r"^(\S+).*\[HSA:(\d+)\]")
    for entry in parse_entries(data_dir / "disease.txt"):
        for line in entry.get("GENE", []):
            match = hsa_pattern.match(line.strip())
            if match:
                mapping.setdefault(match.group(2), match.group(1))

    return mapping


def build_index(df: pd.DataFrame, src_type: str, dst_type: str) -> dict:
    """kegg_relations.csv에서 (src_type -> dst_type) 엣지만 골라 source_id -> {target_id, ...} 인덱스를 만든다."""
    sub = df[(df["source_type"] == src_type) & (df["target_type"] == dst_type)]
    idx = defaultdict(set)
    for s, t in zip(sub["source_id"], sub["target_id"]):
        idx[s].add(t)
    return idx


def load_name_maps(processed_dir: Path) -> dict:
    disease_df = pd.read_csv(processed_dir / "kegg_disease.csv", dtype=str).fillna("")
    drug_df = pd.read_csv(processed_dir / "kegg_drug.csv", dtype=str).fillna("")
    network_df = pd.read_csv(processed_dir / "kegg_network.csv", dtype=str).fillna("")
    dgroup_df = pd.read_csv(processed_dir / "kegg_dgroup.csv", dtype=str).fillna("")

    return {
        "disease_name": dict(zip(disease_df["entry_id"], disease_df["name"])),
        "drug_name": dict(zip(drug_df["entry_id"], drug_df["name"])),
        "drug_dgroup": dict(zip(drug_df["entry_id"], drug_df["dgroup_id"])),
        "network_name": dict(zip(network_df["entry_id"], network_df["name"])),
        "dgroup_name": dict(zip(dgroup_df["entry_id"], dgroup_df["name"])),
        "disease_ids": list(disease_df["entry_id"]),
    }


def short(text: str, n: int = 60) -> str:
    text = text.split(";")[0].strip()
    return text if len(text) <= n else text[: n - 1] + "…"
