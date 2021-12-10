from sarpy.read import readSAR
from sarpy.local_velocity import readEOFParams


tif_path = "C:/Users/张少/Desktop/11.tiff"
eof_path = "C:/Users/张少/Desktop/S1A_OPER_AUX_POEORB_OPOD_20210309T055829_V20151022T225943_20151024T005943.eof"
amplitude, phase = readSAR(tif_path)
all_parameters = readEOFParams(eof_path)

print(amplitude)
print(phase)
print(all_parameters)