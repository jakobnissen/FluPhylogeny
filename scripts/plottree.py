import ete3

pdfpath = snakemake.output[0]
newickpath = snakemake.input[0]
fnapath = snakemake.params[0]

with open(fnapath) as file:
    names = set()
    for line in file:
        if line.startswith('>'):
            names.add(line[1:].strip())

with open(newickpath) as file:
    tree = ete3.Tree(file.read())

style_def = ete3.NodeStyle()
style_obs = ete3.NodeStyle()
style_obs["bgcolor"] = "#ff5555"

for n in tree.traverse():
    if n.is_leaf():
        n.set_style(style_obs if n.name in names else style_def)

tree.render(pdfpath)