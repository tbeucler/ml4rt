"""Contains list of input arguments for training a neural net."""

import os.path
from ml4rt.utils import example_utils
from ml4rt.utils import normalization
from ml4rt.machine_learning import neural_net

TIME_FORMAT = '%Y-%m-%d-%H%M%S'

HOME_DIR_NAME = os.path.expanduser('~')
DEFAULT_EXAMPLE_DIR_NAME = '{0:s}/ml4rt_project/examples/with_new_data'.format(
    HOME_DIR_NAME
)
DEFAULT_NORM_FILE_NAME = '{0:s}/learning_examples_20170101-20181224.nc'.format(
    DEFAULT_EXAMPLE_DIR_NAME
)

DEFAULT_PREDICTOR_NAMES = example_utils.ALL_VECTOR_PREDICTOR_NAMES + [
    example_utils.ZENITH_ANGLE_NAME, example_utils.ALBEDO_NAME,
    example_utils.LATITUDE_NAME, example_utils.LONGITUDE_NAME
]

NET_TYPE_TO_DEFAULT_HEIGHTS_M_AGL = {
    neural_net.U_NET_TYPE_STRING: example_utils.DEFAULT_HEIGHTS_M_AGL[:64],
    neural_net.DENSE_NET_TYPE_STRING: example_utils.DEFAULT_HEIGHTS_M_AGL,
    neural_net.CNN_TYPE_STRING: example_utils.DEFAULT_HEIGHTS_M_AGL
}

NET_TYPE_ARG_NAME = 'net_type_string'
EXAMPLE_DIR_ARG_NAME = 'input_example_dir_name'
INPUT_MODEL_FILE_ARG_NAME = 'input_model_file_name'
OUTPUT_MODEL_DIR_ARG_NAME = 'output_model_dir_name'
PREDICTOR_NAMES_ARG_NAME = 'predictor_names'
TARGET_NAMES_ARG_NAME = 'target_names'
HEIGHTS_ARG_NAME = 'heights_m_agl'
OMIT_HEATING_RATE_ARG_NAME = 'omit_heating_rate'
FIRST_TRAIN_TIME_ARG_NAME = 'first_training_time_string'
LAST_TRAIN_TIME_ARG_NAME = 'last_training_time_string'
FIRST_VALIDN_TIME_ARG_NAME = 'first_validn_time_string'
LAST_VALIDN_TIME_ARG_NAME = 'last_validn_time_string'
NORMALIZATION_FILE_ARG_NAME = 'input_normalization_file_name'
PREDICTOR_NORM_TYPE_ARG_NAME = 'predictor_norm_type_string'
PREDICTOR_MIN_VALUE_ARG_NAME = 'predictor_min_norm_value'
PREDICTOR_MAX_VALUE_ARG_NAME = 'predictor_max_norm_value'
TARGET_NORM_TYPE_ARG_NAME = 'target_norm_type_string'
TARGET_MIN_VALUE_ARG_NAME = 'target_min_norm_value'
TARGET_MAX_VALUE_ARG_NAME = 'target_max_norm_value'
BATCH_SIZE_ARG_NAME = 'num_examples_per_batch'
NUM_EPOCHS_ARG_NAME = 'num_epochs'
NUM_TRAINING_BATCHES_ARG_NAME = 'num_training_batches_per_epoch'
NUM_VALIDN_BATCHES_ARG_NAME = 'num_validn_batches_per_epoch'

NET_TYPE_HELP_STRING = (
    'Neural-net type (must be accepted by `neural_net.check_net_type`).'
)
EXAMPLE_DIR_HELP_STRING = (
    'Name of directory with examples for both training and validation.  Files '
    'therein will be found by `example_io.find_file` and read by '
    '`example_io.read_file`.'
)
INPUT_MODEL_FILE_HELP_STRING = (
    'Path to file with untrained model (defining architecture, optimizer, and '
    'loss function).  Will be read by `neural_net.read_model`.'
)
OUTPUT_MODEL_DIR_HELP_STRING = (
    'Name of output directory.  Model will be saved here.'
)
PREDICTOR_NAMES_HELP_STRING = (
    'List of predictor variables.  Each must be accepted by '
    '`example_utils.check_field_name`.'
)
TARGET_NAMES_HELP_STRING = (
    'List of target variables.  Each must be accepted by '
    '`example_utils.check_field_name`.'
)
HEIGHTS_HELP_STRING = (
    'List of heights (metres above ground level) for profile (vector) '
    'variables.'
)
OMIT_HEATING_RATE_HELP_STRING = (
    'Boolean flag.  If 1, heating rate will not be explicitly predicted by the '
    'model, even if it is in the list of target variables.'
)
TRAIN_TIME_HELP_STRING = (
    'Time (format "yyyy-mm-dd-HHMMSS").  The training period will be '
    '`{0:s}`...`{1:s}`.'
).format(FIRST_TRAIN_TIME_ARG_NAME, LAST_TRAIN_TIME_ARG_NAME)

VALIDN_TIME_HELP_STRING = (
    'Time (format "yyyy-mm-dd-HHMMSS").  The validation (monitoring) period '
    'will be `{0:s}`...`{1:s}`.'
).format(FIRST_VALIDN_TIME_ARG_NAME, LAST_VALIDN_TIME_ARG_NAME)

NORMALIZATION_FILE_HELP_STRING = (
    'Path to normalization file.  Will be read by `example_io.read_file`.'
)
PREDICTOR_NORM_TYPE_HELP_STRING = (
    'Normalization type for predictors (must be accepted by '
    '`normalization._check_normalization_type`).  If you do not want to '
    'normalize, make this an empty string ("").'
)
PREDICTOR_MIN_VALUE_HELP_STRING = (
    'Minimum value if you have chosen min-max normalization.'
)
PREDICTOR_MAX_VALUE_HELP_STRING = (
    'Max value if you have chosen min-max normalization.'
)
TARGET_NORM_TYPE_HELP_STRING = (
    'Same as `{0:s}` but for target variables.'
).format(PREDICTOR_NORM_TYPE_ARG_NAME)

TARGET_MIN_VALUE_HELP_STRING = (
    'Same as `{0:s}` but for target variables.'
).format(PREDICTOR_MIN_VALUE_ARG_NAME)

TARGET_MAX_VALUE_HELP_STRING = (
    'Same as `{0:s}` but for target variables.'
).format(PREDICTOR_MAX_VALUE_ARG_NAME)

BATCH_SIZE_HELP_STRING = (
    'Number of examples in each training and validation (monitoring) batch.'
)
NUM_EPOCHS_HELP_STRING = 'Number of epochs.'
NUM_TRAINING_BATCHES_HELP_STRING = 'Number of training batches per epoch.'
NUM_VALIDN_BATCHES_HELP_STRING = (
    'Number of validation (monitoring) batches per epoch.'
)


def add_input_args(parser_object):
    """Adds input args to ArgumentParser object.

    :param parser_object: Instance of `argparse.ArgumentParser` (may already
        contain some input args).
    :return: parser_object: Same as input but with new args added.
    """

    parser_object.add_argument(
        '--' + NET_TYPE_ARG_NAME, type=str, required=True,
        help=NET_TYPE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + EXAMPLE_DIR_ARG_NAME, type=str, required=False,
        default=DEFAULT_EXAMPLE_DIR_NAME, help=EXAMPLE_DIR_HELP_STRING
    )
    parser_object.add_argument(
        '--' + INPUT_MODEL_FILE_ARG_NAME, type=str, required=True,
        help=INPUT_MODEL_FILE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + OUTPUT_MODEL_DIR_ARG_NAME, type=str, required=True,
        help=OUTPUT_MODEL_DIR_HELP_STRING
    )
    parser_object.add_argument(
        '--' + PREDICTOR_NAMES_ARG_NAME, type=str, nargs='+', required=False,
        default=DEFAULT_PREDICTOR_NAMES, help=PREDICTOR_NAMES_HELP_STRING
    )
    parser_object.add_argument(
        '--' + TARGET_NAMES_ARG_NAME, type=str, nargs='+', required=True,
        help=TARGET_NAMES_HELP_STRING
    )
    parser_object.add_argument(
        '--' + HEIGHTS_ARG_NAME, type=float, nargs='+', required=False,
        default=[-1], help=HEIGHTS_HELP_STRING
    )
    parser_object.add_argument(
        '--' + OMIT_HEATING_RATE_ARG_NAME, type=int, required=False, default=0,
        help=OMIT_HEATING_RATE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + FIRST_TRAIN_TIME_ARG_NAME, type=str, required=False,
        default='2017-01-01-000000', help=TRAIN_TIME_HELP_STRING
    )
    parser_object.add_argument(
        '--' + LAST_TRAIN_TIME_ARG_NAME, type=str, required=False,
        default='2018-12-24-235959', help=TRAIN_TIME_HELP_STRING
    )
    parser_object.add_argument(
        '--' + FIRST_VALIDN_TIME_ARG_NAME, type=str, required=False,
        default='2019-01-01-000000', help=VALIDN_TIME_HELP_STRING
    )
    parser_object.add_argument(
        '--' + LAST_VALIDN_TIME_ARG_NAME, type=str, required=False,
        default='2019-12-24-235959', help=VALIDN_TIME_HELP_STRING
    )
    parser_object.add_argument(
        '--' + NORMALIZATION_FILE_ARG_NAME, type=str, required=False,
        default=DEFAULT_NORM_FILE_NAME, help=NORMALIZATION_FILE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + PREDICTOR_NORM_TYPE_ARG_NAME, type=str, required=False,
        default=normalization.Z_SCORE_NORM_STRING,
        help=PREDICTOR_NORM_TYPE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + PREDICTOR_MIN_VALUE_ARG_NAME, type=float, required=False,
        default=0., help=PREDICTOR_MIN_VALUE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + PREDICTOR_MAX_VALUE_ARG_NAME, type=float, required=False,
        default=1., help=PREDICTOR_MAX_VALUE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + TARGET_NORM_TYPE_ARG_NAME, type=str, required=True,
        help=TARGET_NORM_TYPE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + TARGET_MIN_VALUE_ARG_NAME, type=float, required=False,
        default=0., help=TARGET_MIN_VALUE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + TARGET_MAX_VALUE_ARG_NAME, type=float, required=False,
        default=1., help=TARGET_MAX_VALUE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + BATCH_SIZE_ARG_NAME, type=int, required=False, default=1024,
        help=BATCH_SIZE_HELP_STRING
    )
    parser_object.add_argument(
        '--' + NUM_EPOCHS_ARG_NAME, type=int, required=False, default=1000,
        help=NUM_EPOCHS_HELP_STRING
    )
    parser_object.add_argument(
        '--' + NUM_TRAINING_BATCHES_ARG_NAME, type=int, required=False,
        default=32, help=NUM_TRAINING_BATCHES_HELP_STRING
    )
    parser_object.add_argument(
        '--' + NUM_VALIDN_BATCHES_ARG_NAME, type=int, required=False,
        default=16, help=NUM_VALIDN_BATCHES_HELP_STRING
    )

    return parser_object