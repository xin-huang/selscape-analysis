import svgutils.transform as sg

human = sg.fromfile(snakemake.input.human)
apes  = sg.fromfile(snakemake.input.apes)

human_width  = float(human.width.replace("pt", ""))
human_height = float(human.height.replace("pt", ""))
apes_width   = float(apes.width.replace("pt", ""))

human_root = human.getroot()
apes_root  = apes.getroot()

apes_root.moveto(human_width, 0)

combined = sg.SVGFigure()
combined.set_size((f"{human_width + apes_width}pt", f"{human_height}pt"))
combined.append([human_root, apes_root])
combined.save(snakemake.output.combined)


