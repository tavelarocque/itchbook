#!/usr/bin/env python3
"""Render docs/src/*.md walkthroughs to book-style PDFs in docs/.

Markup understood (deliberately small):
  # Title                 first line
  _subtitle_              second non-blank line, optional
  ## Section / ### Sub    headings
  ```                     fenced code block, rendered verbatim
  | a | b |               pipe table (first row is the header; the |---| row is skipped)
  ![caption](path)        figure, path relative to the repo root, scaled to the text width
  - item                  bullet list item (continuation lines are indented)
  **bold** *italic* `code`
  blank line              paragraph break

Usage: python3 scripts/walkthrough.py [docs/src/FILE.md ...]   (needs reportlab)
"""

from __future__ import annotations

import re
import sys
from html import escape
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Image, KeepTogether, Paragraph, Preformatted, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parent.parent
FONT_DIRS = [Path("/usr/share/fonts/truetype/dejavu"), Path("/Library/Fonts"), Path.home() / "Library/Fonts"]


def font(name: str) -> Path:
    for d in FONT_DIRS:
        if (d / name).exists():
            return d / name
    raise SystemExit(f"font {name} not found; install DejaVu fonts")


for alias, file in [("Serif", "DejaVuSerif.ttf"), ("Serif-Bold", "DejaVuSerif-Bold.ttf"),
                    ("Serif-Italic", "DejaVuSerif-Italic.ttf"), ("Sans-Bold", "DejaVuSans-Bold.ttf"),
                    ("Mono", "DejaVuSansMono.ttf")]:
    pdfmetrics.registerFont(TTFont(alias, font(file)))
pdfmetrics.registerFontFamily("Serif", normal="Serif", bold="Serif-Bold", italic="Serif-Italic")

INK = colors.HexColor("#1b1f24")
MUTED = colors.HexColor("#5b636e")
ACCENT = colors.HexColor("#0f5e9c")
CODE_BG = colors.HexColor("#f3f4f6")
RULE = colors.HexColor("#d0d4da")

S = {
    "title": ParagraphStyle("title", fontName="Sans-Bold", fontSize=21, leading=26, textColor=INK, spaceAfter=4),
    "subtitle": ParagraphStyle("subtitle", fontName="Serif-Italic", fontSize=11, leading=15, textColor=MUTED,
                               spaceAfter=16),
    "h2": ParagraphStyle("h2", fontName="Sans-Bold", fontSize=13.5, leading=17, textColor=ACCENT, spaceBefore=14,
                         spaceAfter=6),
    "h3": ParagraphStyle("h3", fontName="Sans-Bold", fontSize=11, leading=14, textColor=INK, spaceBefore=9,
                         spaceAfter=4),
    "body": ParagraphStyle("body", fontName="Serif", fontSize=10.5, leading=15.4, textColor=INK,
                           alignment=TA_JUSTIFY, spaceAfter=8),
    "bullet": ParagraphStyle("bullet", fontName="Serif", fontSize=10.5, leading=15.4, textColor=INK, leftIndent=16,
                             bulletIndent=5, spaceAfter=3),
    "code": ParagraphStyle("code", fontName="Mono", fontSize=8.2, leading=10.6, textColor=INK),
    "cell": ParagraphStyle("cell", fontName="Serif", fontSize=9, leading=11.5, textColor=INK),
    "cellh": ParagraphStyle("cellh", fontName="Sans-Bold", fontSize=8.8, leading=11.5, textColor=INK),
    "caption": ParagraphStyle("caption", fontName="Serif-Italic", fontSize=9, leading=12, textColor=MUTED,
                              alignment=TA_JUSTIFY, spaceAfter=12),
}


def figure(path: str, caption: str):
    img = Image(str(ROOT / path))
    scale = min(1.0, 6.5 * inch / img.imageWidth)
    img.drawWidth, img.drawHeight = img.imageWidth * scale, img.imageHeight * scale
    img._restrictSize(6.5 * inch, 4.6 * inch)
    return [KeepTogether([img, Spacer(1, 4), Paragraph(inline(caption), S["caption"])])]


def inline(text: str) -> str:
    parts = re.split(r"(`[^`]+`)", text)
    out = []
    for p in parts:
        if p.startswith("`") and p.endswith("`") and len(p) > 1:
            out.append(f'<font face="Mono" size="9">{escape(p[1:-1])}</font>')
        else:
            p = escape(p)
            p = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", p)
            p = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<i>\1</i>", p)
            out.append(p)
    return "".join(out)


def code_block(lines: list[str]):
    t = Table([[Preformatted("\n".join(lines), S["code"])]], colWidths=[6.5 * inch])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), CODE_BG), ("LEFTPADDING", (0, 0), (-1, -1), 8),
                           ("RIGHTPADDING", (0, 0), (-1, -1), 8), ("TOPPADDING", (0, 0), (-1, -1), 6),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
    return [t, Spacer(1, 8)]


def table(rows: list[str]):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows
             if not re.fullmatch(r"\|?[\s:|-]+\|?", r.strip())]
    data = [[Paragraph(inline(c), S["cellh" if i == 0 else "cell"]) for c in row] for i, row in enumerate(cells)]
    n = max(len(r) for r in data)
    t = Table(data, colWidths=[6.5 * inch / n] * n, repeatRows=1)
    t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, RULE),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"), ("TOPPADDING", (0, 0), (-1, -1), 3),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
    return [t, Spacer(1, 10)]


def render(src: Path) -> Path:
    lines = src.read_text().splitlines()
    story, para, i = [], [], 0

    def flush():
        if para:
            story.append(Paragraph(inline(" ".join(para)), S["body"]))
            para.clear()

    title = lines[0].removeprefix("# ").strip()
    story.append(Paragraph(inline(title), S["title"]))
    i = 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i < len(lines) and lines[i].startswith("_") and lines[i].rstrip().endswith("_"):
        story.append(Paragraph(inline(lines[i].strip()[1:-1]), S["subtitle"]))
        i += 1
    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            flush()
            j = i + 1
            while j < len(lines) and not lines[j].startswith("```"):
                j += 1
            story += code_block(lines[i + 1:j])
            i = j + 1
            continue
        fig = re.fullmatch(r"!\[(.*)\]\((.+)\)", line.strip())
        if fig:
            flush()
            story += figure(fig.group(2), fig.group(1))
            i += 1
            continue
        if line.startswith("|"):
            flush()
            j = i
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            story += table(lines[i:j])
            i = j
            continue
        if line.startswith("- "):
            flush()
            item = [line[2:].strip()]
            while i + 1 < len(lines) and lines[i + 1].startswith("  ") and lines[i + 1].strip():
                i += 1
                item.append(lines[i].strip())
            story.append(Paragraph(inline(" ".join(item)), S["bullet"], bulletText="\u2022"))
            i += 1
            if i >= len(lines) or not lines[i].startswith("- "):
                story.append(Spacer(1, 5))
            continue
        if line.startswith("### "):
            flush()
            story.append(Paragraph(inline(line[4:]), S["h3"]))
        elif line.startswith("## "):
            flush()
            story.append(Paragraph(inline(line[3:]), S["h2"]))
        elif not line.strip():
            flush()
        else:
            para.append(line.strip())
        i += 1
    flush()

    out = ROOT / "docs" / (src.stem + ".pdf")

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Serif-Italic", 8.5)
        canvas.setFillColor(MUTED)
        canvas.drawString(inch, 0.6 * inch, f"itchbook · {title}")
        canvas.drawRightString(letter[0] - inch, 0.6 * inch, str(doc.page))
        canvas.restoreState()

    doc = SimpleDocTemplate(str(out), pagesize=letter, leftMargin=inch, rightMargin=inch, topMargin=0.9 * inch,
                            bottomMargin=0.9 * inch, title=title, author="Tave La Rocque")
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return out


if __name__ == "__main__":
    srcs = [Path(a) for a in sys.argv[1:]] or sorted((ROOT / "docs" / "src").glob("*.md"))
    for s in srcs:
        print(render(s))
