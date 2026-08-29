#!/usr/bin/env python
"""intake/network 스케치의 전사를 재현한다 (L0 감사용).

무엇을 하는가:
  ① EXIF orientation 보정 + 90° CCW 회전 (원본은 눕혀 찍힘)
  ② 그림에서 센 비드 21개의 위치·반지름을 **박아넣고** 원본 위에 겹쳐 그린다
  ③ 접촉 판정으로 링크를 재구성해 위상(고리 수·분기점·자유단)을 출력한다

★ 비드 좌표는 자동 검출이 아니라 **사람(=이 스크립트를 쓴 세션)이 눈으로 읽은 값**이다.
  자동 Hough 검출은 이 환경에 opencv/skimage 가 없어 못 썼다. 그래서 이 스크립트의
  역할은 "세었다"가 아니라 **"센 결과를 남에게 검증받을 수 있게 만든다"** 이다 —
  out/network_sketch_annotated.png 를 원본과 나란히 보면 빠뜨림·중복이 눈에 보인다.

  전사 규칙(skill bd-intake §1): 입자 수가 적으면(≲50) 정확히 세고, 많으면 세지 말고
  적절히 많은 N 을 제안한다. 여기서는 21개라 **정확히 셌다**. 다만 그 21이
  시뮬레이션의 N 이라는 뜻은 아니다 — observation.yaml A6 참조.

실행:
  $PY scratch/read_network_sketch.py
"""
from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "intake/network/KakaoTalk_Photo_2026-08-06-14-36-42.jpeg"
OUT = ROOT / "intake/network/transcription_annotated.png"

# 회전·크롭 규약: exif_transpose → ROTATE_90(CCW) → crop(260,600,1920,2260)
# 이 좌표계(1660x1660)에서 읽은 비드 (x, y, r).  r 은 눈대중 — 접촉 판정에만 쓴다.
CROP = (260, 600, 1920, 2260)
BEADS = [
    (700, 260, 150), (795, 510, 155),                                   # 위쪽 가지
    (1070, 665, 95), (1215, 570, 80), (1335, 520, 85),                  # 우상 가지
    (940, 730, 105), (1015, 975, 105), (775, 865, 75), (650, 975, 80),  # 중앙
    (505, 1075, 65), (425, 1035, 45), (370, 1015, 35),                  # 좌 가지
    (305, 975, 55), (185, 950, 72),
    (395, 1205, 62), (320, 1330, 58),                                   # 좌하 가지
    (1105, 1175, 100), (1225, 1150, 65), (1345, 1055, 62), (1465, 1000, 68),  # 우하 가지
    (1160, 1420, 100),                                                  # 아래 가지
]
# 링크(접한 원 쌍)도 **눈으로 읽은 데이터**다. 1-based 번호쌍, 위 BEADS 순서 기준.
# ⚠️ 자동(중심거리 < r_i+r_j) 재구성은 이 좌표로 **작동하지 않는다** — 아래 참조.
LINKS_READ = [
    (1, 2), (2, 6), (6, 3), (3, 4), (4, 5),          # 위쪽 + 우상 가지
    (6, 7), (6, 8),                                   # 6 = 4갈래 분기점
    (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14),   # 좌 가지
    (10, 15), (15, 16),                               # 10 = 3갈래 분기점
    (7, 17), (17, 18), (18, 19), (19, 20), (17, 21),  # 17 = 3갈래 분기점
]


def contact_ratio(beads, i, j):
    """중심거리 / (r_i + r_j). 1 미만이면 원이 겹친다."""
    xi, yi, ri = beads[i]
    xj, yj, rj = beads[j]
    return ((xi - xj) ** 2 + (yi - yj) ** 2) ** 0.5 / (ri + rj)


def main() -> int:
    if not SRC.exists():
        print(f"원본 없음: {SRC}", file=sys.stderr)
        return 1

    img = ImageOps.exif_transpose(Image.open(SRC)).transpose(Image.ROTATE_90)
    dia = img.crop(CROP).convert("RGB")

    lk = [(a - 1, b - 1) for a, b in LINKS_READ]
    deg = [0] * len(BEADS)
    for i, j in lk:
        deg[i] += 1
        deg[j] += 1

    # 연결성분 (union-find 없이 BFS 한 번 — 21개짜리라 충분)
    adj = {i: [] for i in range(len(BEADS))}
    for i, j in lk:
        adj[i].append(j)
        adj[j].append(i)
    seen, comps = set(), 0
    for s in range(len(BEADS)):
        if s in seen:
            continue
        comps += 1
        stack = [s]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(v for v in adj[u] if v not in seen)

    n, e = len(BEADS), len(lk)
    loops = e - n + comps          # 순환 차원 (Betti-1)

    print(f"비드            {n}")
    print(f"링크            {e}   (그림에서 읽음)")
    print(f"연결성분        {comps}")
    print(f"독립 고리       {loops}   ← 0 이면 트리 (응력 나르는 닫힌 경로 없음)")
    print(f"자유단(deg=1)   {sum(1 for d in deg if d == 1)}")
    print(f"3갈래(deg=3)    {sum(1 for d in deg if d == 3)}")
    print(f"4갈래(deg=4)    {sum(1 for d in deg if d == 4)}")
    print(f"평균 배위수 z   {2 * e / n:.3f}")

    # ── 자동 재구성은 왜 못 쓰는가 (TOL 을 맞추는 대신 **재서 보고한다**) ──────
    lset = {frozenset(p) for p in lk}
    lr = sorted(contact_ratio(BEADS, i, j) for i, j in lk)
    nr = sorted(contact_ratio(BEADS, i, j)
                for i in range(n) for j in range(i + 1, n)
                if frozenset((i, j)) not in lset)
    print("\n중심거리/(r_i+r_j) — 눈대중 반지름의 품질 진단")
    print(f"  링크인 쌍  ({len(lr):2d}개)  최소 {lr[0]:.2f} · 중앙 {lr[len(lr)//2]:.2f} · 최대 {lr[-1]:.2f}")
    print(f"  아닌 쌍   ({len(nr):3d}개)  최소 {nr[0]:.2f}")
    print(f"  → 링크 최대({lr[-1]:.2f}) {'>' if lr[-1] > nr[0] else '<'} 비링크 최소({nr[0]:.2f}) "
          f": 거리 하나로는 {'분리 불가' if lr[-1] > nr[0] else '분리 가능'}")
    print("  ★ BEADS 의 r 은 **번호를 얹으려고 작게** 잡은 값이라 그려진 타원의 반축이")
    print("    아니다. 그래서 접촉 판정을 자동화하면 링크를 놓친다 — TOL 을 올려")
    print("    맞추면 원하는 답에 끼워맞추는 것이 되므로, 링크는 눈으로 읽고")
    print("    이 진단만 남긴다 (CLAUDE.md 규칙 6).")

    dr = ImageDraw.Draw(dia)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 52)
    except OSError:
        font = None
    for i, j in lk:                                    # 링크 먼저 (비드 아래로)
        dr.line([BEADS[i][:2], BEADS[j][:2]], fill=(0, 120, 255), width=7)
    for idx, (x, y, r) in enumerate(BEADS):
        dr.ellipse([x - r, y - r, x + r, y + r], outline=(255, 0, 0), width=5)
        dr.text((x - 18, y - 26), str(idx + 1), fill=(0, 150, 0), font=font)
    dia.save(OUT)
    print(f"\n주석 그림 → {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
