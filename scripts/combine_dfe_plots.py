# Copyright 2025 Xin Huang and Simon Chen
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


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


