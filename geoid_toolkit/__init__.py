import geoid_toolkit.compute
import geoid_toolkit.datum
import geoid_toolkit.interpolate
import geoid_toolkit.math
import geoid_toolkit.spatial
import geoid_toolkit.utilities
import geoid_toolkit.version
from geoid_toolkit.read_ICGEM_harmonics import read_ICGEM_harmonics
from geoid_toolkit.read_topography_harmonics import read_topography_harmonics

# import functions for backwards compatibility
from geoid_toolkit.compute import (
    geoid_undulation,
    gravity_anomaly,
    gravity_disturbance,
    height_anomaly,
    real_potential,
    topographic_potential,
)
from geoid_toolkit.datum import (
    ref_ellipsoid,
    norm_gravity,
    norm_potential,
)

# executable scripts
from geoid_toolkit import scripts

# get version number
__version__ = geoid_toolkit.version.version
