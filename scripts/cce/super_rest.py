import scri
from scri.SpEC.file_io import write_to_h5
import numpy as np

# time after junk radiation
t0 = 1000
# roughly 2 initial orbital periods
w0 = 800
# interpolate to these times
t_interp = np.linspace(0, 2500, 16384)

# Read the data
abd = scri.create_abd_from_h5(
    file_name="CharacteristicExtractReduction.h5",
    file_format="spectrecce_v1",
    ch_mass=1.0,
    t_0_superrest=t0,
    padding_time=w0,
    t_interpolate=t_interp
)

# Write data
write_to_h5(abd.h, "Strain.h5")
