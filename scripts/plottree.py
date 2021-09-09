import ete3

segment = snakemake.wildcards["segment"]
flutype = snakemake.wildcards["flutype"]
pdfpath = snakemake.output[0]
newickpath = snakemake.input[0]

with open(f"subtypes/{segment}.txt") as file:
    names = set()
    lines = map(str.split, filter(None, map(str.strip, file)))
    for (sample, seqname, _flutype, _) in lines:
        if flutype == _flutype:
            names.add(seqname)

with open(newickpath) as file:
    tree = ete3.Tree(file.read())

style_def = ete3.NodeStyle()
style_obs = ete3.NodeStyle()
style_obs["bgcolor"] = "#ff5555"

for n in tree.traverse():
    if n.is_leaf():
        n.set_style(style_obs if n.name in names else style_def)

tree.render(pdfpath)