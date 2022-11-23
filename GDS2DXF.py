import pya
import sys

def GDS2DXF(GDSfile, DXFfile):
    layout = pya.Layout()
    layout.read(GDSfile)
    layout.write(DXFfile)

def KeepOnlyTop(GDSfile):
    layout = pya.Layout()
    layout.read(GDSfile)
    cells = layout.each_cell()

    cellIndexs = []
    for cell in cells:
        if cell.basic_name() != 'top':
            cellIndexs = cellIndexs + [cell.cell_index()]
    layout.delete_cells(cellIndexs)
    layout.write(GDSfile)

def FaltTop(GDSfile):
    layout = pya.Layout()
    layout.read(GDSfile)
    layout.flatten(layout.cell('top').cell_index(), 100, 1)
    layout.write(GDSfile)


if __name__ == '__main__':
    FaltTop(sys.argv[1])
    KeepOnlyTop(sys.argv[1])
    GDS2DXF(sys.argv[1], sys.argv[2])
