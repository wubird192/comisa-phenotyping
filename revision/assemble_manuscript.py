#!/usr/bin/env python3
"""assemble_manuscript.py — build the full revised manuscript (Markdown) from ms_parts/ and the PLOS reference list.

Usage: python3 -I assemble_manuscript.py <ms_parts_dir> <References_PLOS_vN.md> <out.md>
"""
import sys, os, re
parts_dir, refs_md, out = sys.argv[1:4]
refs = open(refs_md, encoding="utf-8").read()
sec1 = refs.split("## 1. ")[1].split("## 2. ")[0]
keymap = dict((m.group(2), int(m.group(1))) for m in re.finditer(r"^(\d+)\. .*?\n   <!-- \[(.+?)\] -->", sec1, re.M))
reflist = re.sub(r"\n   <!-- \[.+?\] -->", "", sec1.split("\n", 2)[2]).strip()
order = ["01_front_abstract_intro.md","02_methods.md","03_results.md","04_discussion.md","05_backmatter.md"]
body = "\n\n".join(open(os.path.join(parts_dir, f), encoding="utf-8").read().strip() for f in order)
body += "\n\n**REFERENCES**\n\n" + reflist + "\n\n"
body += open(os.path.join(parts_dir, "06_tables.md"), encoding="utf-8").read().strip() + "\n\n"
legends = open(os.path.join(parts_dir, "07_figure_legends.md"), encoding="utf-8").read()
for n in (3, 2, 1):
    m = re.search(r"^\*\*Fig %d\..*$" % n, legends, re.M)
    if not m: continue
    paras = body.split("\n\n")
    for i, para in enumerate(paras):
        if re.search(r"\(Fig %d\)|\bFig %d\b" % (n, n), para) and not para.startswith("**Fig"):
            paras.insert(i + 1, m.group(0)); break
    body = "\n\n".join(paras)
body += "\n"
unknown, used = [], set()
def repl(m):
    keys = [k.strip() for k in m.group(1).split(";")]
    nums = []
    for k in keys:
        if k in keymap: nums.append(keymap[k]); used.add(k)
        elif k.startswith("cited "): return m.group(0)
        else: unknown.append(k); nums.append(f"?{k}")
    return "[" + ",".join(str(n) for n in nums) + "]"
body = re.sub(r"\[((?:[A-Za-z][^\[\]]*?\d{4}[^\[\]]*?)(?:;\s*[^\[\]]+?)*)\]", repl, body)
open(out, "w", encoding="utf-8").write(body)
print("wrote", out, len(body.split()), "words")
print("unknown keys:", sorted(set(unknown)) or "none")
print("unused references:", sorted(set(keymap) - used, key=lambda k: keymap[k]) or "none")
