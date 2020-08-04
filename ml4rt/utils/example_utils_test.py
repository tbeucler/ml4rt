"""Unit tests for example_utils.py."""

import copy
import unittest
import numpy
from ml4rt.utils import example_utils

TOLERANCE = 1e-6

# The following constants are used to test _find_nonzero_runs.
FIRST_VALUES = numpy.full(10, 0.)
FIRST_START_INDICES = numpy.array([], dtype=int)
FIRST_END_INDICES = numpy.array([], dtype=int)

SECOND_VALUES = numpy.random.uniform(low=-5., high=5., size=10)
SECOND_VALUES[numpy.absolute(SECOND_VALUES) < TOLERANCE] = 1.
SECOND_START_INDICES = numpy.array([0], dtype=int)
SECOND_END_INDICES = numpy.array([9], dtype=int)

THIRD_VALUES = numpy.array([
    0, 0, 0, 0.2, 1.2, 3.5, 13.4, 31.2, 56.8, 90.1, 129.0, 172.3, 219.1, 249.4,
    263.7, 0, 0, 0
])
THIRD_START_INDICES = numpy.array([3], dtype=int)
THIRD_END_INDICES = numpy.array([14], dtype=int)

FOURTH_VALUES = numpy.concatenate((THIRD_VALUES, THIRD_VALUES, THIRD_VALUES))
FOURTH_START_INDICES = numpy.array([3, 21, 39], dtype=int)
FOURTH_END_INDICES = numpy.array([14, 32, 50], dtype=int)

# The following constants are used to test get_grid_cell_edges and
# get_grid_cell_widths.
CENTER_HEIGHTS_M_AGL = numpy.array([
    10, 20, 40, 60, 80, 100, 30000, 33000, 36000, 39000, 42000, 46000, 50000
], dtype=float)

EDGE_HEIGHTS_M_AGL = numpy.array([
    5, 15, 30, 50, 70, 90, 15050, 31500, 34500, 37500, 40500, 44000, 48000,
    52000
], dtype=float)

GRID_CELL_WIDTHS_METRES = numpy.array([
    10, 15, 20, 20, 20, 14960, 16450, 3000, 3000, 3000, 3500, 4000, 4000
], dtype=float)

# The following constants are used to test fluxes_to_heating_rate,
# fluxes_actual_to_increments, and fluxes_increments_to_actual.
THIS_UP_FLUX_MATRIX_W_M02 = numpy.array([
    [100, 150, 200, 250, 300, 350],
    [400, 500, 600, 700, 800, 900],
    [0, 0, 0, 0, 0, 0]
], dtype=float)

THIS_DOWN_FLUX_MATRIX_W_M02 = numpy.array([
    [50, 125, 200, 275, 350, 425],
    [500, 550, 600, 650, 700, 750],
    [1000, 1000, 1000, 1000, 1000, 1000]
], dtype=float)

THIS_PRESSURE_MATRIX_PASCALS = 100 * numpy.array([
    [1000, 950, 900, 850, 800, 750],
    [1000, 900, 800, 700, 600, 500],
    [1000, 950, 900, 850, 800, 750]
], dtype=float)

THIS_PRESSURE_DIFF_MATRIX_PASCALS = 100 * numpy.array([
    [-50, -50, -50, -50, -50, -50],
    [-100, -100, -100, -100, -100, -100],
    [-50, -50, -50, -50, -50, -50]
], dtype=float)

# THIS_NET_FLUX_MATRIX_W_M02 = numpy.array([
#     [-50, -25, 0, 25, 50, 75],
#     [100, 50, 0, -50, -100, -150],
#     [1000, 1000, 1000, 1000, 1000, 1000]
# ], dtype=float)

THIS_NET_FLUX_DIFF_MATRIX_W02 = numpy.array([
    [25, 25, 25, 25, 25, 25],
    [-50, -50, -50, -50, -50, -50],
    [0, 0, 0, 0, 0, 0]
], dtype=float)

THESE_HEIGHTS_M_AGL = numpy.array([25, 50, 100, 500, 1000, 5000], dtype=float)
VALID_TIMES_UNIX_SEC = numpy.array([300, 600, 900], dtype=int)

THIS_WIDTH_MATRIX_METRES = numpy.array([
    [25, 37.5, 225, 450, 2250, 4000],
    [25, 37.5, 225, 450, 2250, 4000],
    [25, 37.5, 225, 450, 2250, 4000]
])

THIS_COEFF = example_utils.DAYS_TO_SECONDS * (
    example_utils.GRAVITY_CONSTANT_M_S02 /
    example_utils.DRY_AIR_SPECIFIC_HEAT_J_KG01_K01
)

THIS_HEATING_RATE_MATRIX_K_DAY01 = THIS_COEFF * (
    THIS_NET_FLUX_DIFF_MATRIX_W02 /
    numpy.absolute(THIS_PRESSURE_DIFF_MATRIX_PASCALS)
)

THIS_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    THIS_PRESSURE_MATRIX_PASCALS, axis=-1
)
THESE_VECTOR_PREDICTOR_NAMES = [example_utils.PRESSURE_NAME]

THIS_VECTOR_TARGET_MATRIX = numpy.stack(
    (THIS_UP_FLUX_MATRIX_W_M02, THIS_DOWN_FLUX_MATRIX_W_M02), axis=-1
)
THESE_VECTOR_TARGET_NAMES = [
    example_utils.SHORTWAVE_UP_FLUX_NAME, example_utils.SHORTWAVE_DOWN_FLUX_NAME
]

EXAMPLE_DICT_FLUXES_ONLY = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX + 0.,
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX + 0.,
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: THESE_HEIGHTS_M_AGL
}

THIS_VECTOR_TARGET_MATRIX = numpy.stack((
    THIS_UP_FLUX_MATRIX_W_M02, THIS_DOWN_FLUX_MATRIX_W_M02,
    THIS_HEATING_RATE_MATRIX_K_DAY01
), axis=-1)

THESE_VECTOR_TARGET_NAMES = [
    example_utils.SHORTWAVE_UP_FLUX_NAME,
    example_utils.SHORTWAVE_DOWN_FLUX_NAME,
    example_utils.SHORTWAVE_HEATING_RATE_NAME
]

EXAMPLE_DICT_WITH_HEATING_RATE = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX + 0.,
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX + 0.,
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: THESE_HEIGHTS_M_AGL
}

THIS_UP_FLUX_INC_MATRIX_W_M02 = numpy.array([
    [100, 50, 50, 50, 50, 50],
    [400, 100, 100, 100, 100, 100],
    [0, 0, 0, 0, 0, 0]
], dtype=float)

THIS_DOWN_FLUX_INC_MATRIX_W_M02 = numpy.array([
    [50, 75, 75, 75, 75, 75],
    [500, 50, 50, 50, 50, 50],
    [1000, 0, 0, 0, 0, 0]
], dtype=float)

THIS_UP_FLUX_INC_MATRIX_W_M03 = (
    THIS_UP_FLUX_INC_MATRIX_W_M02 / THIS_WIDTH_MATRIX_METRES
)
THIS_DOWN_FLUX_INC_MATRIX_W_M03 = (
    THIS_DOWN_FLUX_INC_MATRIX_W_M02 / THIS_WIDTH_MATRIX_METRES
)

THIS_VECTOR_TARGET_MATRIX = numpy.stack((
    THIS_UP_FLUX_MATRIX_W_M02, THIS_DOWN_FLUX_MATRIX_W_M02,
    THIS_DOWN_FLUX_INC_MATRIX_W_M03, THIS_UP_FLUX_INC_MATRIX_W_M03
), axis=-1)

THESE_VECTOR_TARGET_NAMES = [
    example_utils.SHORTWAVE_UP_FLUX_NAME,
    example_utils.SHORTWAVE_DOWN_FLUX_NAME,
    example_utils.SHORTWAVE_DOWN_FLUX_INC_NAME,
    example_utils.SHORTWAVE_UP_FLUX_INC_NAME
]

EXAMPLE_DICT_WITH_INCREMENTS = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX + 0.,
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX + 0.,
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: THESE_HEIGHTS_M_AGL
}

# The following constants are used to test find_cloud_layers.
THIS_LWP_MATRIX_KG_M02 = 0.001 * numpy.array([
    [0, 0, 0, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6, 6, 6, 6],
    [10, 20, 50, 70, 90, 90, 90, 90, 90, 90, 90, 110, 110, 130, 130, 150, 150,
     210],
    [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51]
], dtype=float)

MIN_PATH_KG_M02 = 0.05

CLOUD_MASK_MATRIX = numpy.array([
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
], dtype=bool)

CLOUD_LAYER_COUNTS = numpy.array([0, 2, 1], dtype=int)

THESE_HEIGHTS_M_AGL = numpy.array([
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
], dtype=float)

THESE_VECTOR_PREDICTOR_NAMES = [example_utils.UPWARD_LIQUID_WATER_PATH_NAME]
THIS_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    THIS_LWP_MATRIX_KG_M02, axis=-1
)

EXAMPLE_DICT_FOR_CLOUD_MASK = {
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: THESE_VECTOR_PREDICTOR_NAMES,
    example_utils.HEIGHTS_KEY: THESE_HEIGHTS_M_AGL
}

# The following constants are used to test concat_examples.
FIRST_TIMES_UNIX_SEC = numpy.array([0, 300, 600, 1200], dtype=int)
FIRST_STANDARD_ATMO_FLAGS = numpy.array([0, 1, 2, 3], dtype=int)
FIRST_EXAMPLE_ID_STRINGS = ['foo', 'bar', 'moo', 'hal']

SCALAR_PREDICTOR_NAMES = [
    example_utils.ZENITH_ANGLE_NAME, example_utils.LATITUDE_NAME
]

FIRST_ZENITH_ANGLES_RADIANS = numpy.array([0, 1, 2, 3], dtype=float)
FIRST_LATITUDES_DEG_N = numpy.array([40.02, 40.02, 40.02, 40.02])
FIRST_SCALAR_PREDICTOR_MATRIX = numpy.transpose(numpy.vstack(
    (FIRST_ZENITH_ANGLES_RADIANS, FIRST_LATITUDES_DEG_N)
))

VECTOR_PREDICTOR_NAMES = [example_utils.TEMPERATURE_NAME]
HEIGHTS_M_AGL = numpy.array([100, 500], dtype=float)

FIRST_TEMP_MATRIX_KELVINS = numpy.array([
    [290, 295],
    [289, 294],
    [288, 293],
    [287, 292.5]
])
FIRST_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    FIRST_TEMP_MATRIX_KELVINS, axis=-1
)

SCALAR_TARGET_NAMES = [example_utils.SHORTWAVE_SURFACE_DOWN_FLUX_NAME]

FIRST_SURFACE_DOWN_FLUXES_W_M02 = numpy.array(
    [200, 200, 200, 200], dtype=float
)
FIRST_SCALAR_TARGET_MATRIX = numpy.reshape(
    FIRST_SURFACE_DOWN_FLUXES_W_M02,
    (len(FIRST_SURFACE_DOWN_FLUXES_W_M02), 1)
)

VECTOR_TARGET_NAMES = [
    example_utils.SHORTWAVE_DOWN_FLUX_NAME, example_utils.SHORTWAVE_UP_FLUX_NAME
]

FIRST_DOWN_FLUX_MATRIX_W_M02 = numpy.array([
    [300, 200],
    [500, 300],
    [450, 450],
    [200, 100]
], dtype=float)

FIRST_UP_FLUX_MATRIX_W_M02 = numpy.array([
    [150, 150],
    [200, 150],
    [300, 350],
    [400, 100]
], dtype=float)

FIRST_VECTOR_TARGET_MATRIX = numpy.stack(
    (FIRST_DOWN_FLUX_MATRIX_W_M02, FIRST_UP_FLUX_MATRIX_W_M02), axis=-1
)

FIRST_EXAMPLE_DICT = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY: FIRST_SCALAR_PREDICTOR_MATRIX,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY: FIRST_VECTOR_PREDICTOR_MATRIX,
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX,
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: FIRST_VECTOR_TARGET_MATRIX,
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC,
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS,
    example_utils.EXAMPLE_IDS_KEY: FIRST_EXAMPLE_ID_STRINGS
}

SECOND_EXAMPLE_ID_STRINGS = ['FOO', 'BAR', 'MOO', 'HAL']

SECOND_EXAMPLE_DICT = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY: FIRST_SCALAR_PREDICTOR_MATRIX * 2,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY: FIRST_VECTOR_PREDICTOR_MATRIX * 3,
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX * 4,
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: FIRST_VECTOR_TARGET_MATRIX * 5,
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC * 6,
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS + 1,
    example_utils.EXAMPLE_IDS_KEY: SECOND_EXAMPLE_ID_STRINGS
}

CONCAT_EXAMPLE_DICT = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY: numpy.concatenate(
        (FIRST_SCALAR_PREDICTOR_MATRIX, FIRST_SCALAR_PREDICTOR_MATRIX * 2),
        axis=0
    ),
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY: numpy.concatenate(
        (FIRST_VECTOR_PREDICTOR_MATRIX, FIRST_VECTOR_PREDICTOR_MATRIX * 3),
        axis=0
    ),
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: numpy.concatenate(
        (FIRST_SCALAR_TARGET_MATRIX, FIRST_SCALAR_TARGET_MATRIX * 4),
        axis=0
    ),
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: numpy.concatenate(
        (FIRST_VECTOR_TARGET_MATRIX, FIRST_VECTOR_TARGET_MATRIX * 5),
        axis=0
    ),
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: numpy.concatenate(
        (FIRST_TIMES_UNIX_SEC, FIRST_TIMES_UNIX_SEC * 6),
        axis=0
    ),
    example_utils.STANDARD_ATMO_FLAGS_KEY: numpy.concatenate(
        (FIRST_STANDARD_ATMO_FLAGS, FIRST_STANDARD_ATMO_FLAGS + 1),
        axis=0
    ),
    example_utils.EXAMPLE_IDS_KEY:
        FIRST_EXAMPLE_ID_STRINGS + SECOND_EXAMPLE_ID_STRINGS
}

# The following constants are used to test subset_by_time.
FIRST_SUBSET_TIME_UNIX_SEC = 1
LAST_SUBSET_TIME_UNIX_SEC = 600
GOOD_INDICES_SELECT_TIMES = numpy.array([1, 2], dtype=int)

FIRST_EXAMPLE_DICT_SELECT_TIMES = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY:
        FIRST_SCALAR_PREDICTOR_MATRIX[1:3, ...],
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        FIRST_VECTOR_PREDICTOR_MATRIX[1:3, ...],
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX[1:3, ...],
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: FIRST_VECTOR_TARGET_MATRIX[1:3, ...],
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC[1:3, ...],
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS[1:3, ...],
    example_utils.EXAMPLE_IDS_KEY: FIRST_EXAMPLE_ID_STRINGS[1:3]
}

# The following constants are used to test subset_by_standard_atmo.
STANDARD_ATMO_ENUM = 2
GOOD_INDICES_STANDARD_ATMO = numpy.array([2], dtype=int)

FIRST_EXAMPLE_DICT_SELECT_ATMO_TYPES = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY:
        FIRST_SCALAR_PREDICTOR_MATRIX[[2], ...],
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        FIRST_VECTOR_PREDICTOR_MATRIX[[2], ...],
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX[[2], ...],
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: FIRST_VECTOR_TARGET_MATRIX[[2], ...],
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC[[2], ...],
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS[[2], ...],
    example_utils.EXAMPLE_IDS_KEY: [FIRST_EXAMPLE_ID_STRINGS[2]]
}

# The following constants are used to test subset_by_field.
FIELD_NAMES_TO_KEEP = [
    example_utils.SHORTWAVE_UP_FLUX_NAME, example_utils.LATITUDE_NAME,
    example_utils.SHORTWAVE_DOWN_FLUX_NAME
]

FIRST_EXAMPLE_DICT_SELECT_FIELDS = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: [SCALAR_PREDICTOR_NAMES[1]],
    example_utils.SCALAR_PREDICTOR_VALS_KEY:
        FIRST_SCALAR_PREDICTOR_MATRIX[..., [1]],
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: [],
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        FIRST_VECTOR_PREDICTOR_MATRIX[..., []],
    example_utils.SCALAR_TARGET_NAMES_KEY: [],
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX[..., []],
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES[::-1],
    example_utils.VECTOR_TARGET_VALS_KEY: FIRST_VECTOR_TARGET_MATRIX[..., ::-1],
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC,
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS,
    example_utils.EXAMPLE_IDS_KEY: FIRST_EXAMPLE_ID_STRINGS
}

# The following constants are used to test subset_by_height.
HEIGHTS_TO_KEEP_M_AGL = numpy.array([500, 100], dtype=float)

FIRST_EXAMPLE_DICT_SELECT_HEIGHTS = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY: FIRST_SCALAR_PREDICTOR_MATRIX,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        FIRST_VECTOR_PREDICTOR_MATRIX[:, ::-1, :],
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: FIRST_SCALAR_TARGET_MATRIX,
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY:
        FIRST_VECTOR_TARGET_MATRIX[:, ::-1, :],
    example_utils.HEIGHTS_KEY: HEIGHTS_TO_KEEP_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC,
    example_utils.STANDARD_ATMO_FLAGS_KEY: FIRST_STANDARD_ATMO_FLAGS,
    example_utils.EXAMPLE_IDS_KEY: FIRST_EXAMPLE_ID_STRINGS
}

# The following constants are used to test find_examples.
ALL_ID_STRINGS = ['south_boulder', 'bear', 'green', 'flagstaff', 'sanitas']
DESIRED_ID_STRINGS_0MISSING = ['green', 'bear']
RELEVANT_INDICES_0MISSING = numpy.array([2, 1], dtype=int)

DESIRED_ID_STRINGS_2MISSING = ['green', 'paiute', 'bear', 'audubon']
RELEVANT_INDICES_2MISSING = numpy.array([2, -1, 1, -1], dtype=int)

# The following constants are used to test subset_by_index.
FIRST_EXAMPLE_DICT_SELECT_INDICES = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY:
        FIRST_SCALAR_PREDICTOR_MATRIX[[2, 1], ...],
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        FIRST_VECTOR_PREDICTOR_MATRIX[[2, 1], ...],
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY:
        FIRST_SCALAR_TARGET_MATRIX[[2, 1], ...],
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY:
        FIRST_VECTOR_TARGET_MATRIX[[2, 1], ...],
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL,
    example_utils.VALID_TIMES_KEY: FIRST_TIMES_UNIX_SEC[[2, 1], ...],
    example_utils.STANDARD_ATMO_FLAGS_KEY:
        FIRST_STANDARD_ATMO_FLAGS[[2, 1], ...],
    example_utils.EXAMPLE_IDS_KEY:
        [FIRST_EXAMPLE_ID_STRINGS[2], FIRST_EXAMPLE_ID_STRINGS[1]]
}

# The following constants are used to test average_examples.
THIS_SCALAR_PREDICTOR_MATRIX = numpy.array([[1.5, 40.02]])
THIS_VECTOR_PREDICTOR_MATRIX = numpy.array([[288.5, 293.625]])
THIS_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    THIS_VECTOR_PREDICTOR_MATRIX, axis=-1
)
THIS_SCALAR_TARGET_MATRIX = numpy.array([[200]], dtype=float)
THIS_VECTOR_TARGET_MATRIX = numpy.array([
    [362.5, 262.5],
    [262.5, 187.5]
])

FIRST_EXAMPLE_DICT_AVERAGE = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: SCALAR_PREDICTOR_NAMES,
    example_utils.SCALAR_PREDICTOR_VALS_KEY: THIS_SCALAR_PREDICTOR_MATRIX,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: VECTOR_PREDICTOR_NAMES,
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX,
    example_utils.SCALAR_TARGET_NAMES_KEY: SCALAR_TARGET_NAMES,
    example_utils.SCALAR_TARGET_VALS_KEY: THIS_SCALAR_TARGET_MATRIX,
    example_utils.VECTOR_TARGET_NAMES_KEY: VECTOR_TARGET_NAMES,
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX,
    example_utils.HEIGHTS_KEY: HEIGHTS_M_AGL
}

# The following constants are used to test create_example_ids and
# parse_example_ids.
LATITUDES_FOR_IDS_DEG_N = numpy.array([40, 40.04, 53.5, 40.0381113])
LONGITUDES_FOR_IDS_DEG_E = numpy.array([255, 254.74, 246.5, 254.7440276])
ZENITH_ANGLES_FOR_IDS_RAD = numpy.array([0.5, 0.666, 0.7777777, 1])
TEMPERATURES_FOR_IDS_KELVINS = numpy.array([230, 240, 250, 260], dtype=float)

TIMES_FOR_IDS_UNIX_SEC = numpy.array([
    0, int(1e7), int(1e8), int(1e9)
], dtype=int)

STANDARD_ATMO_FLAGS_FOR_IDS = numpy.array([
    example_utils.MIDLATITUDE_WINTER_ENUM,
    example_utils.MIDLATITUDE_WINTER_ENUM,
    example_utils.SUBARCTIC_WINTER_ENUM,
    example_utils.MIDLATITUDE_WINTER_ENUM
], dtype=int)

THIS_SCALAR_PREDICTOR_MATRIX = numpy.transpose(numpy.vstack((
    LATITUDES_FOR_IDS_DEG_N, LONGITUDES_FOR_IDS_DEG_E, ZENITH_ANGLES_FOR_IDS_RAD
)))

THIS_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    TEMPERATURES_FOR_IDS_KELVINS, axis=-1
)
THIS_VECTOR_PREDICTOR_MATRIX = numpy.expand_dims(
    THIS_VECTOR_PREDICTOR_MATRIX, axis=-1
)

EXAMPLE_DICT_FOR_IDS = {
    example_utils.SCALAR_PREDICTOR_NAMES_KEY: [
        example_utils.LATITUDE_NAME, example_utils.LONGITUDE_NAME,
        example_utils.ZENITH_ANGLE_NAME
    ],
    example_utils.SCALAR_PREDICTOR_VALS_KEY: THIS_SCALAR_PREDICTOR_MATRIX,
    example_utils.VECTOR_PREDICTOR_NAMES_KEY: [example_utils.TEMPERATURE_NAME],
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX,
    example_utils.HEIGHTS_KEY: numpy.array([10.]),
    example_utils.VALID_TIMES_KEY: TIMES_FOR_IDS_UNIX_SEC,
    example_utils.STANDARD_ATMO_FLAGS_KEY: STANDARD_ATMO_FLAGS_FOR_IDS
}

EXAMPLE_ID_STRINGS = [
    'lat=40.000000_long=255.000000_zenith-angle-rad=0.500000_'
    'time=0000000000_atmo={0:d}_temp-10m-kelvins=230.000000'.format(
        example_utils.MIDLATITUDE_WINTER_ENUM
    ),
    'lat=40.040000_long=254.740000_zenith-angle-rad=0.666000_'
    'time=0010000000_atmo={0:d}_temp-10m-kelvins=240.000000'.format(
        example_utils.MIDLATITUDE_WINTER_ENUM
    ),
    'lat=53.500000_long=246.500000_zenith-angle-rad=0.777778_'
    'time=0100000000_atmo={0:d}_temp-10m-kelvins=250.000000'.format(
        example_utils.SUBARCTIC_WINTER_ENUM
    ),
    'lat=40.038111_long=254.744028_zenith-angle-rad=1.000000_'
    'time=1000000000_atmo={0:d}_temp-10m-kelvins=260.000000'.format(
        example_utils.MIDLATITUDE_WINTER_ENUM
    )
]

# The following constants are used to test create_fake_heights.
REAL_HEIGHTS_M_AGL = numpy.array(
    [10, 20, 40, 60, 80, 10000, 50000], dtype=float
)
NUM_PADDING_HEIGHTS = 10
PADDED_HEIGHTS_M_AGL = numpy.array([
    10, 20, 40, 60, 80, 10000, 50000, 1050000, 2050000, 3050000, 4050000,
    5050000, 6050000, 7050000, 8050000, 9050000, 10050000
], dtype=float)

# The following constants are used to test _add_height_padding.
THIS_HUMIDITY_MATRIX_KG_KG01 = 0.001 * numpy.array([
    [1, 2, 3, 4, 5, 6, 7],
    [2, 4, 6, 8, 10, 12, 14],
    [3, 6, 9, 12, 15, 18, 21]
], dtype=float)

THIS_TEMPERATURE_MATRIX_KELVINS = 273.15 + numpy.array([
    [10, 11, 12, 13, 14, 15, 16],
    [20, 21, 22, 23, 24, 25, 26],
    [30, 31, 32, 33, 34, 35, 36]
], dtype=float)

THESE_VECTOR_PREDICTOR_NAMES = [
    example_utils.SPECIFIC_HUMIDITY_NAME, example_utils.TEMPERATURE_NAME
]
THIS_VECTOR_PREDICTOR_MATRIX = numpy.stack((
    THIS_HUMIDITY_MATRIX_KG_KG01, THIS_TEMPERATURE_MATRIX_KELVINS,
), axis=-1)

THIS_UP_FLUX_MATRIX_W_M02 = numpy.array([
    [100, 150, 200, 250, 300, 350, 400],
    [400, 500, 600, 700, 800, 900, 1000],
    [0, 0, 0, 0, 0, 0, 0]
], dtype=float)

THIS_DOWN_FLUX_MATRIX_W_M02 = numpy.array([
    [50, 125, 200, 275, 350, 425, 525],
    [500, 550, 600, 650, 700, 750, 850],
    [1000, 1000, 1000, 1000, 1000, 1000, 1400]
], dtype=float)

THESE_VECTOR_TARGET_NAMES = [
    example_utils.SHORTWAVE_UP_FLUX_NAME, example_utils.SHORTWAVE_DOWN_FLUX_NAME
]
THIS_VECTOR_TARGET_MATRIX = numpy.stack(
    (THIS_UP_FLUX_MATRIX_W_M02, THIS_DOWN_FLUX_MATRIX_W_M02), axis=-1
)

EXAMPLE_DICT_SANS_PADDING = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX + 0.,
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX + 0.,
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: REAL_HEIGHTS_M_AGL + 0.
}

THIS_HUMIDITY_MATRIX_KG_KG01 = 0.001 * numpy.array([
    [1, 2, 3, 4, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [2, 4, 6, 8, 10, 12, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14],
    [3, 6, 9, 12, 15, 18, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21]
], dtype=float)

THIS_TEMPERATURE_MATRIX_KELVINS = 273.15 + numpy.array([
    [10, 11, 12, 13, 14, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16],
    [20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26],
    [30, 31, 32, 33, 34, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36]
], dtype=float)

THIS_VECTOR_PREDICTOR_MATRIX = numpy.stack((
    THIS_HUMIDITY_MATRIX_KG_KG01, THIS_TEMPERATURE_MATRIX_KELVINS,
), axis=-1)

THIS_UP_FLUX_MATRIX_W_M02 = numpy.array([
    [100, 150, 200, 250, 300, 350, 400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [400, 500, 600, 700, 800, 900, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
], dtype=float)

THIS_DOWN_FLUX_MATRIX_W_M02 = numpy.array([
    [50, 125, 200, 275, 350, 425, 525, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [500, 550, 600, 650, 700, 750, 850, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1000, 1000, 1000, 1000, 1000, 1000, 1400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
], dtype=float)

THIS_VECTOR_TARGET_MATRIX = numpy.stack(
    (THIS_UP_FLUX_MATRIX_W_M02, THIS_DOWN_FLUX_MATRIX_W_M02), axis=-1
)

EXAMPLE_DICT_WITH_PADDING = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY: THIS_VECTOR_PREDICTOR_MATRIX + 0.,
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX + 0.,
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: PADDED_HEIGHTS_M_AGL + 0.
}

# The following constants are used to test subset_by_height with padding.
PADDED_HEIGHTS_TO_KEEP_M_AGL = numpy.array([
    2050000, 3050000, 4050000, 5050000, 6050000, 7050000, 8050000, 9050000,
    10050000
], dtype=float)

EXAMPLE_DICT_PADDED_SELECT_HEIGHTS = {
    example_utils.VECTOR_PREDICTOR_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_PREDICTOR_NAMES),
    example_utils.VECTOR_PREDICTOR_VALS_KEY:
        THIS_VECTOR_PREDICTOR_MATRIX[:, -9:, :],
    example_utils.VECTOR_TARGET_NAMES_KEY:
        copy.deepcopy(THESE_VECTOR_TARGET_NAMES),
    example_utils.VECTOR_TARGET_VALS_KEY: THIS_VECTOR_TARGET_MATRIX[:, -9:, :],
    example_utils.VALID_TIMES_KEY: VALID_TIMES_UNIX_SEC,
    example_utils.HEIGHTS_KEY: PADDED_HEIGHTS_TO_KEEP_M_AGL + 0.
}


def _compare_example_dicts(first_example_dict, second_example_dict):
    """Compares two dictionaries with learning examples.

    :param first_example_dict: See doc for `example_io.read_file`.
    :param second_example_dict: Same.
    :return: are_dicts_equal: Boolean flag.
    """

    first_keys = list(first_example_dict.keys())
    second_keys = list(first_example_dict.keys())
    if set(first_keys) != set(second_keys):
        return False

    keys_to_compare = [
        example_utils.SCALAR_PREDICTOR_VALS_KEY,
        example_utils.VECTOR_PREDICTOR_VALS_KEY,
        example_utils.SCALAR_TARGET_VALS_KEY,
        example_utils.VECTOR_TARGET_VALS_KEY,
        example_utils.HEIGHTS_KEY
    ]

    for this_key in keys_to_compare:
        if this_key not in first_example_dict:
            continue

        if not numpy.allclose(
                first_example_dict[this_key], second_example_dict[this_key],
                atol=TOLERANCE
        ):
            return False

    try:
        if not numpy.array_equal(
                first_example_dict[example_utils.VALID_TIMES_KEY],
                second_example_dict[example_utils.VALID_TIMES_KEY]
        ):
            return False
    except KeyError:
        pass

    keys_to_compare = [
        example_utils.SCALAR_PREDICTOR_NAMES_KEY,
        example_utils.VECTOR_PREDICTOR_NAMES_KEY,
        example_utils.SCALAR_TARGET_NAMES_KEY,
        example_utils.VECTOR_TARGET_NAMES_KEY,
        example_utils.EXAMPLE_IDS_KEY
    ]

    for this_key in keys_to_compare:
        if this_key not in first_example_dict:
            continue

        if first_example_dict[this_key] != second_example_dict[this_key]:
            return False

    return True


class ExampleUtilsTests(unittest.TestCase):
    """Each method is a unit test for example_utils.py."""

    def test_find_nonzero_runs_first(self):
        """Ensures correct output from find_nonzero_runs.

        In this case, using first set of values.
        """

        these_start_indices, these_end_indices = (
            example_utils._find_nonzero_runs(FIRST_VALUES)
        )
        self.assertTrue(numpy.array_equal(
            these_start_indices, FIRST_START_INDICES
        ))
        self.assertTrue(numpy.array_equal(
            these_end_indices, FIRST_END_INDICES
        ))

    def test_find_nonzero_runs_second(self):
        """Ensures correct output from find_nonzero_runs.

        In this case, using second set of values.
        """

        these_start_indices, these_end_indices = (
            example_utils._find_nonzero_runs(SECOND_VALUES)
        )
        self.assertTrue(numpy.array_equal(
            these_start_indices, SECOND_START_INDICES
        ))
        self.assertTrue(numpy.array_equal(
            these_end_indices, SECOND_END_INDICES
        ))

    def test_find_nonzero_runs_third(self):
        """Ensures correct output from find_nonzero_runs.

        In this case, using third set of values.
        """

        these_start_indices, these_end_indices = (
            example_utils._find_nonzero_runs(THIRD_VALUES)
        )
        self.assertTrue(numpy.array_equal(
            these_start_indices, THIRD_START_INDICES
        ))
        self.assertTrue(numpy.array_equal(
            these_end_indices, THIRD_END_INDICES
        ))

    def test_find_nonzero_runs_fourth(self):
        """Ensures correct output from find_nonzero_runs.

        In this case, using fourth set of values.
        """

        these_start_indices, these_end_indices = (
            example_utils._find_nonzero_runs(FOURTH_VALUES)
        )
        self.assertTrue(numpy.array_equal(
            these_start_indices, FOURTH_START_INDICES
        ))
        self.assertTrue(numpy.array_equal(
            these_end_indices, FOURTH_END_INDICES
        ))

    def test_get_grid_cell_edges(self):
        """Ensures correct output from get_grid_cell_edges."""

        these_edge_heights_m_agl = (
            example_utils.get_grid_cell_edges(CENTER_HEIGHTS_M_AGL)
        )
        self.assertTrue(numpy.allclose(
            these_edge_heights_m_agl, EDGE_HEIGHTS_M_AGL, atol=TOLERANCE
        ))

    def test_get_grid_cell_widths(self):
        """Ensures correct output from get_grid_cell_widths."""

        these_widths_metres = (
            example_utils.get_grid_cell_widths(EDGE_HEIGHTS_M_AGL)
        )
        self.assertTrue(numpy.allclose(
            these_widths_metres, GRID_CELL_WIDTHS_METRES, atol=TOLERANCE
        ))

    def test_fluxes_to_heating_rate(self):
        """Ensures correct output from fluxes_to_heating_rate."""

        this_example_dict = example_utils.fluxes_to_heating_rate(
            copy.deepcopy(EXAMPLE_DICT_FLUXES_ONLY)
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, EXAMPLE_DICT_WITH_HEATING_RATE
        ))

    def test_fluxes_actual_to_increments(self):
        """Ensures correct output from fluxes_actual_to_increments."""

        this_example_dict = example_utils.fluxes_actual_to_increments(
            copy.deepcopy(EXAMPLE_DICT_FLUXES_ONLY)
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, EXAMPLE_DICT_WITH_INCREMENTS
        ))

    def test_fluxes_increments_to_actual(self):
        """Ensures correct output from fluxes_increments_to_actual."""

        this_example_dict = example_utils.fluxes_increments_to_actual(
            copy.deepcopy(EXAMPLE_DICT_WITH_INCREMENTS)
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, EXAMPLE_DICT_WITH_INCREMENTS
        ))

    def test_find_cloud_layers(self):
        """Ensures correct output from find_cloud_layers."""

        this_mask_matrix, these_cloud_layer_counts = (
            example_utils.find_cloud_layers(
                example_dict=EXAMPLE_DICT_FOR_CLOUD_MASK,
                min_path_kg_m02=MIN_PATH_KG_M02, for_ice=False
            )
        )

        self.assertTrue(numpy.array_equal(this_mask_matrix, CLOUD_MASK_MATRIX))
        self.assertTrue(numpy.array_equal(
            these_cloud_layer_counts, CLOUD_LAYER_COUNTS
        ))

    def test_concat_examples_good(self):
        """Ensures correct output from concat_examples.

        In this case, not expecting an error.
        """

        this_example_dict = example_utils.concat_examples(
            [FIRST_EXAMPLE_DICT, SECOND_EXAMPLE_DICT]
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, CONCAT_EXAMPLE_DICT
        ))

    def test_concat_examples_bad_heights(self):
        """Ensures correct output from concat_examples.

        In this case, expecting an error due to mismatched heights.
        """

        this_second_example_dict = copy.deepcopy(SECOND_EXAMPLE_DICT)
        this_second_example_dict[example_utils.HEIGHTS_KEY] += 1

        with self.assertRaises(ValueError):
            example_utils.concat_examples(
                [FIRST_EXAMPLE_DICT, this_second_example_dict]
            )

    def test_concat_examples_bad_fields(self):
        """Ensures correct output from concat_examples.

        In this case, expecting an error due to mismatched fields.
        """

        this_second_example_dict = copy.deepcopy(SECOND_EXAMPLE_DICT)

        this_second_example_dict[
            example_utils.SCALAR_PREDICTOR_NAMES_KEY
        ].append(example_utils.ALBEDO_NAME)

        with self.assertRaises(ValueError):
            example_utils.concat_examples(
                [FIRST_EXAMPLE_DICT, this_second_example_dict]
            )

    def test_get_field_zenith_no_height(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for zenith angle at no particular height.
        """

        this_vector = example_utils.get_field_from_dict(
            example_dict=FIRST_EXAMPLE_DICT,
            field_name=example_utils.ZENITH_ANGLE_NAME, height_m_agl=None
        )

        self.assertTrue(numpy.allclose(
            this_vector, FIRST_ZENITH_ANGLES_RADIANS
        ))

    def test_get_field_zenith_with_height(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for zenith angle at particular height.
        """

        this_vector = example_utils.get_field_from_dict(
            example_dict=FIRST_EXAMPLE_DICT,
            field_name=example_utils.ZENITH_ANGLE_NAME, height_m_agl=10.
        )

        self.assertTrue(numpy.allclose(
            this_vector, FIRST_ZENITH_ANGLES_RADIANS
        ))

    def test_get_field_temperature_no_height(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for temperature at no particular height.
        """

        this_matrix = example_utils.get_field_from_dict(
            example_dict=FIRST_EXAMPLE_DICT,
            field_name=example_utils.TEMPERATURE_NAME, height_m_agl=None
        )

        self.assertTrue(numpy.allclose(this_matrix, FIRST_TEMP_MATRIX_KELVINS))

    def test_get_field_temperature_100m(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for temperature at 100 m AGL.
        """

        this_vector = example_utils.get_field_from_dict(
            example_dict=FIRST_EXAMPLE_DICT,
            field_name=example_utils.TEMPERATURE_NAME, height_m_agl=100.
        )

        self.assertTrue(numpy.allclose(
            this_vector, FIRST_TEMP_MATRIX_KELVINS[:, 0]
        ))

    def test_get_field_temperature_500m(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for temperature at 500 m AGL.
        """

        this_vector = example_utils.get_field_from_dict(
            example_dict=FIRST_EXAMPLE_DICT,
            field_name=example_utils.TEMPERATURE_NAME, height_m_agl=500.
        )

        self.assertTrue(numpy.allclose(
            this_vector, FIRST_TEMP_MATRIX_KELVINS[:, 1]
        ))

    def test_get_field_temperature_600m(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for temperature at 600 m AGL (unavailable).
        """

        with self.assertRaises(ValueError):
            example_utils.get_field_from_dict(
                example_dict=FIRST_EXAMPLE_DICT,
                field_name=example_utils.TEMPERATURE_NAME, height_m_agl=600.
            )

    def test_get_field_lwp(self):
        """Ensures correct output from get_field_from_dict.

        In this case, looking for liquid-water path (unavailable).
        """

        with self.assertRaises(ValueError):
            example_utils.get_field_from_dict(
                example_dict=FIRST_EXAMPLE_DICT,
                field_name=example_utils.LIQUID_WATER_PATH_NAME,
                height_m_agl=None
            )

    def test_subset_by_time(self):
        """Ensures correct output from subset_by_time."""

        this_example_dict, these_indices = example_utils.subset_by_time(
            example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT),
            first_time_unix_sec=FIRST_SUBSET_TIME_UNIX_SEC,
            last_time_unix_sec=LAST_SUBSET_TIME_UNIX_SEC
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_SELECT_TIMES
        ))
        self.assertTrue(numpy.array_equal(
            these_indices, GOOD_INDICES_SELECT_TIMES
        ))

    def test_subset_by_standard_atmo(self):
        """Ensures correct output from subset_by_standard_atmo."""

        this_example_dict, these_indices = (
            example_utils.subset_by_standard_atmo(
                example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT),
                standard_atmo_enum=STANDARD_ATMO_ENUM
            )
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_SELECT_ATMO_TYPES
        ))
        self.assertTrue(numpy.array_equal(
            these_indices, GOOD_INDICES_STANDARD_ATMO
        ))

    def test_subset_by_field(self):
        """Ensures correct output from subset_by_field."""

        this_example_dict = example_utils.subset_by_field(
            example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT),
            field_names=FIELD_NAMES_TO_KEEP
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_SELECT_FIELDS
        ))

    def test_subset_by_height(self):
        """Ensures correct output from subset_by_height."""

        this_example_dict = example_utils.subset_by_height(
            example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT),
            heights_m_agl=HEIGHTS_TO_KEEP_M_AGL
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_SELECT_HEIGHTS
        ))

    def test_find_examples_0missing(self):
        """Ensures correct output from find_examples.

        In this case, no desired examples are missing.
        """

        these_indices = example_utils.find_examples(
            all_id_strings=ALL_ID_STRINGS,
            desired_id_strings=DESIRED_ID_STRINGS_0MISSING, allow_missing=False
        )
        self.assertTrue(numpy.array_equal(
            these_indices, RELEVANT_INDICES_0MISSING
        ))

    def test_find_examples_2missing_allowed(self):
        """Ensures correct output from find_examples.

        In this case, 2 desired examples are missing but this is allowed.
        """

        these_indices = example_utils.find_examples(
            all_id_strings=ALL_ID_STRINGS,
            desired_id_strings=DESIRED_ID_STRINGS_2MISSING, allow_missing=True
        )
        self.assertTrue(numpy.array_equal(
            these_indices, RELEVANT_INDICES_2MISSING
        ))

    def test_find_examples_2missing_disallowed(self):
        """Ensures correct output from find_examples.

        In this case, 2 desired examples are missing and this is *not* allowed.
        """

        with self.assertRaises(ValueError):
            example_utils.find_examples(
                all_id_strings=ALL_ID_STRINGS,
                desired_id_strings=DESIRED_ID_STRINGS_2MISSING,
                allow_missing=False
            )

    def test_subset_by_index(self):
        """Ensures correct output from subset_by_index."""

        this_example_dict = example_utils.subset_by_index(
            example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT),
            desired_indices=RELEVANT_INDICES_0MISSING
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_SELECT_INDICES
        ))

    def test_average_examples(self):
        """Ensures correct output from average_examples."""

        this_example_dict = example_utils.average_examples(
            example_dict=copy.deepcopy(FIRST_EXAMPLE_DICT), use_pmm=False
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, FIRST_EXAMPLE_DICT_AVERAGE
        ))

    def test_create_example_ids(self):
        """Ensures correct output from create_example_ids."""

        these_id_strings = example_utils.create_example_ids(
            EXAMPLE_DICT_FOR_IDS
        )
        self.assertTrue(these_id_strings == EXAMPLE_ID_STRINGS)

    def test_parse_example_ids(self):
        """Ensures correct output from parse_example_ids."""

        metadata_dict = example_utils.parse_example_ids(EXAMPLE_ID_STRINGS)
        these_latitudes_deg_n = metadata_dict[example_utils.LATITUDES_KEY]
        these_longitudes_deg_e = metadata_dict[example_utils.LONGITUDES_KEY]
        these_zenith_angles_rad = metadata_dict[example_utils.ZENITH_ANGLES_KEY]
        these_times_unix_sec = metadata_dict[example_utils.VALID_TIMES_KEY]
        these_standard_atmo_flags = (
            metadata_dict[example_utils.STANDARD_ATMO_FLAGS_KEY]
        )
        these_10m_temps_kelvins = (
            metadata_dict[example_utils.TEMPERATURES_10M_KEY]
        )

        self.assertTrue(numpy.allclose(
            these_latitudes_deg_n, LATITUDES_FOR_IDS_DEG_N, atol=TOLERANCE
        ))
        self.assertTrue(numpy.allclose(
            these_longitudes_deg_e, LONGITUDES_FOR_IDS_DEG_E, atol=TOLERANCE
        ))
        self.assertTrue(numpy.allclose(
            these_zenith_angles_rad, ZENITH_ANGLES_FOR_IDS_RAD, atol=TOLERANCE
        ))
        self.assertTrue(numpy.array_equal(
            these_times_unix_sec, TIMES_FOR_IDS_UNIX_SEC
        ))
        self.assertTrue(numpy.array_equal(
            these_standard_atmo_flags, STANDARD_ATMO_FLAGS_FOR_IDS
        ))
        self.assertTrue(numpy.allclose(
            these_10m_temps_kelvins, TEMPERATURES_FOR_IDS_KELVINS,
            atol=TOLERANCE
        ))

    def test_create_fake_heights(self):
        """Ensures correct output from create_fake_heights."""

        these_heights_m_agl = example_utils.create_fake_heights(
            real_heights_m_agl=REAL_HEIGHTS_M_AGL,
            num_padding_heights=NUM_PADDING_HEIGHTS
        )

        self.assertTrue(numpy.allclose(
            these_heights_m_agl, PADDED_HEIGHTS_M_AGL, atol=TOLERANCE
        ))

    def test_add_height_padding(self):
        """Ensures correct output from _add_height_padding."""

        this_example_dict = example_utils._add_height_padding(
            example_dict=copy.deepcopy(EXAMPLE_DICT_SANS_PADDING),
            desired_heights_m_agl=PADDED_HEIGHTS_M_AGL
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, EXAMPLE_DICT_WITH_PADDING
        ))

    def test_subset_by_height_with_padding(self):
        """Ensures correct output from subset_by_height.

        In this case, some of the desired heights are not there yet (and should
        be added at the top of the profile).
        """

        this_example_dict = example_utils.subset_by_height(
            example_dict=copy.deepcopy(EXAMPLE_DICT_SANS_PADDING),
            heights_m_agl=PADDED_HEIGHTS_TO_KEEP_M_AGL
        )

        self.assertTrue(_compare_example_dicts(
            this_example_dict, EXAMPLE_DICT_PADDED_SELECT_HEIGHTS
        ))


if __name__ == '__main__':
    unittest.main()
