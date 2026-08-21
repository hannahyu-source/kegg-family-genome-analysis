"""KEGG flat file(ENTRY ... /// 형식) 공용 파서.

KEGG flat file은 필드명이 앞 12칸에 오고 내용이 13번째 칸부터 시작하는
고정폭 포맷을 쓴다. 필드명이 비어 있는 줄(앞 12칸이 공백)은 바로 위
필드의 연속 줄로 취급한다. 엔트리는 "///" 줄로 구분된다.
"""

import re
from pathlib import Path

FIELD_WIDTH = 12


def parse_entries(path: Path):
    """flat file을 읽어 엔트리 리스트를 반환한다.

    각 엔트리는 {field_name: [content_line, ...]} 형태의 dict.
    필드 순서는 파일에 등장한 순서 그대로 보존된다(dict 삽입 순서).
    """
    entries = []
    current = None
    last_field = None

    with open(path, encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")

            if line.strip() == "///":
                if current:
                    entries.append(current)
                current = None
                last_field = None
                continue

            if line.strip() == "":
                continue

            if current is None:
                current = {}

            field_part = line[:FIELD_WIDTH]
            content_part = line[FIELD_WIDTH:]
            field_name = field_part.strip()

            if field_name:
                last_field = field_name
                current.setdefault(field_name, []).append(content_part)
            elif last_field:
                current.setdefault(last_field, []).append(content_part)

    if current:
        entries.append(current)

    return entries


def entry_id_and_type(entry: dict):
    """ENTRY 필드에서 (id, type)을 뽑는다. 예: 'DG00001   DGroup' -> ('DG00001', 'DGroup')."""
    lines = entry.get("ENTRY", [])
    if not lines:
        return None, None
    parts = re.split(r"\s{2,}", lines[0].strip(), maxsplit=1)
    entry_id = parts[0] if parts else None
    entry_type = parts[1].strip() if len(parts) > 1 else None
    return entry_id, entry_type


def field_text(entry: dict, field: str, sep: str = " ") -> str:
    """필드의 모든 줄을 하나의 문자열로 합친다."""
    return sep.join(line.strip() for line in entry.get(field, []) if line.strip())


def field_lines(entry: dict, field: str):
    """필드의 원본 줄 리스트(좌우 공백 제거)를 반환한다."""
    return [line.strip() for line in entry.get(field, []) if line.strip()]
