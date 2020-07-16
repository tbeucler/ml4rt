"""Plots saliency maps (one for each example)."""

import os.path
import argparse
import numpy
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot
from gewittergefahr.gg_utils import file_system_utils
from ml4rt.io import prediction_io
from ml4rt.io import example_io
from ml4rt.machine_learning import saliency
from ml4rt.machine_learning import neural_net
from ml4rt.plotting import profile_plotting

SEPARATOR_STRING = '\n\n' + '*' * 50 + '\n\n'

METRES_TO_KM = 0.001

ORANGE_COLOUR = numpy.array([217, 95, 2], dtype=float) / 255
PURPLE_COLOUR = numpy.array([117, 112, 179], dtype=float) / 255
GREEN_COLOUR = numpy.array([27, 158, 119], dtype=float) / 255

FIRST_VECTOR_PREDICTOR_NAMES = [
    example_io.TEMPERATURE_NAME, example_io.SPECIFIC_HUMIDITY_NAME,
    example_io.PRESSURE_NAME
]
FIRST_VECTOR_PREDICTOR_COLOURS = [ORANGE_COLOUR, PURPLE_COLOUR, GREEN_COLOUR]

SECOND_VECTOR_PREDICTOR_NAMES = [
    example_io.LIQUID_WATER_CONTENT_NAME, example_io.LIQUID_WATER_PATH_NAME,
    example_io.UPWARD_LIQUID_WATER_PATH_NAME
]
SECOND_VECTOR_PREDICTOR_COLOURS = [ORANGE_COLOUR, PURPLE_COLOUR, GREEN_COLOUR]

THIRD_VECTOR_PREDICTOR_NAMES = [
    example_io.ICE_WATER_CONTENT_NAME, example_io.ICE_WATER_PATH_NAME,
    example_io.UPWARD_ICE_WATER_PATH_NAME
]
THIRD_VECTOR_PREDICTOR_COLOURS = [ORANGE_COLOUR, PURPLE_COLOUR, GREEN_COLOUR]

FIRST_SCALAR_PREDICTOR_NAMES = [
    example_io.ZENITH_ANGLE_NAME, example_io.ALBEDO_NAME,
    example_io.COLUMN_LIQUID_WATER_PATH_NAME
]
FIRST_SCALAR_PREDICTOR_COLOURS = [ORANGE_COLOUR, PURPLE_COLOUR, GREEN_COLOUR]

SECOND_SCALAR_PREDICTOR_NAMES = [
    example_io.LATITUDE_NAME, example_io.LONGITUDE_NAME,
    example_io.COLUMN_ICE_WATER_PATH_NAME
]
SECOND_SCALAR_PREDICTOR_COLOURS = [ORANGE_COLOUR, PURPLE_COLOUR, GREEN_COLOUR]

PREDICTOR_NAME_TO_VERBOSE = {
    example_io.TEMPERATURE_NAME: 'Temperature',
    example_io.SPECIFIC_HUMIDITY_NAME: 'Specific\nhumidity',
    example_io.PRESSURE_NAME: 'Pressure',
    example_io.LIQUID_WATER_CONTENT_NAME: 'LWC',
    example_io.ICE_WATER_CONTENT_NAME: 'IWC',
    example_io.LIQUID_WATER_PATH_NAME: 'Downward LWP',
    example_io.ICE_WATER_PATH_NAME: 'Downward IWP',
    example_io.UPWARD_LIQUID_WATER_PATH_NAME: 'Upward LWP',
    example_io.UPWARD_ICE_WATER_PATH_NAME: 'Upward IWP',
    example_io.ZENITH_ANGLE_NAME: 'Zenith angle',
    example_io.LATITUDE_NAME: 'Latitude',
    example_io.LONGITUDE_NAME: 'Longitude',
    example_io.ALBEDO_NAME: 'Albedo',
    example_io.COLUMN_LIQUID_WATER_PATH_NAME: 'Column LWP',
    example_io.COLUMN_ICE_WATER_PATH_NAME: 'Column IWP'
}

TARGET_NAME_TO_VERBOSE = {
    example_io.SHORTWAVE_HEATING_RATE_NAME: 'heating rate',
    example_io.SHORTWAVE_UP_FLUX_NAME: 'up flux',
    example_io.SHORTWAVE_DOWN_FLUX_NAME: 'down flux',
    example_io.SHORTWAVE_TOA_UP_FLUX_NAME: 'TOA up flux',
    example_io.SHORTWAVE_SURFACE_DOWN_FLUX_NAME: 'sfc down flux'
}

TARGET_NAME_TO_UNITS = {
    example_io.SHORTWAVE_HEATING_RATE_NAME: r'K day$^{-1}$',
    example_io.SHORTWAVE_UP_FLUX_NAME: r'W m$^{-2}$',
    example_io.SHORTWAVE_DOWN_FLUX_NAME: r'W m$^{-2}$',
    example_io.SHORTWAVE_TOA_UP_FLUX_NAME: r'W m$^{-2}$',
    example_io.SHORTWAVE_SURFACE_DOWN_FLUX_NAME: r'W m$^{-2}$'
}

LEGEND_BOUNDING_BOX_DICT = {
    'facecolor': 'white',
    'alpha': 0.7,
    'edgecolor': 'black',
    'linewidth': 2,
    'boxstyle': 'round'
}

REFERENCE_LINE_COLOUR = numpy.full(3, 152. / 255)
FIGURE_WIDTH_INCHES = 15
FIGURE_HEIGHT_INCHES = 15
FIGURE_RESOLUTION_DPI = 300

DEFAULT_FONT_SIZE = 24
LEGEND_FONT_SIZE = 16

pyplot.rc('font', size=DEFAULT_FONT_SIZE)
pyplot.rc('axes', titlesize=DEFAULT_FONT_SIZE)
pyplot.rc('axes', labelsize=DEFAULT_FONT_SIZE)
pyplot.rc('xtick', labelsize=DEFAULT_FONT_SIZE)
pyplot.rc('ytick', labelsize=DEFAULT_FONT_SIZE)
pyplot.rc('legend', fontsize=DEFAULT_FONT_SIZE)
pyplot.rc('figure', titlesize=DEFAULT_FONT_SIZE)

SALIENCY_FILE_ARG_NAME = 'input_saliency_file_name'
PREDICTION_FILE_ARG_NAME = 'input_prediction_file_name'
OUTPUT_DIR_ARG_NAME = 'output_dir_name'

SALIENCY_FILE_HELP_STRING = (
    'Path to saliency file (will be read by `saliency.read_standard_file`).'
)
PREDICTION_FILE_HELP_STRING = (
    'Path to prediction file (will be read by `prediction_io.read_file`).  For '
    'each example in the saliency file, this script will find the predicted/'
    'actual target values in the prediction file and include these in the '
    'legend of the saliency plot.  If saliency file contains values for '
    'non-output neuron, this file is not needed.'
)
OUTPUT_DIR_HELP_STRING = (
    'Name of output directory (figures will be saved here).'
)

INPUT_ARG_PARSER = argparse.ArgumentParser()
INPUT_ARG_PARSER.add_argument(
    '--' + SALIENCY_FILE_ARG_NAME, type=str, required=True,
    help=SALIENCY_FILE_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + PREDICTION_FILE_ARG_NAME, type=str, required=False, default='',
    help=PREDICTION_FILE_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + OUTPUT_DIR_ARG_NAME, type=str, required=True,
    help=OUTPUT_DIR_HELP_STRING
)


def _get_target_values(
        prediction_file_name, model_metadata_dict, example_id_strings,
        target_field_name, target_height_m_agl):
    """Returns predicted and actual target values.

    E = number of examples

    :param prediction_file_name: See documentation at top of file.
    :param model_metadata_dict: Dictionary returned by
        `neural_net.read_metafile`.
    :param example_id_strings: length-E list of example IDs.  Will return target
        values only for these examples.
    :param target_field_name: Name of target variable.
    :param target_height_m_agl: Height of target variable (metres above ground
        level).
    :return: predicted_values: length-E numpy array of predicted target values.
    :return: actual_values: length-E numpy array of actual target values.
    """

    print((
        'Reading predicted and actual target values from: "{0:s}"...'
    ).format(
        prediction_file_name
    ))
    prediction_dict = prediction_io.read_file(prediction_file_name)

    example_indices = numpy.array([
        prediction_dict[prediction_io.EXAMPLE_IDS_KEY].index(id)
        for id in example_id_strings
    ], dtype=int)

    generator_option_dict = (
        model_metadata_dict[neural_net.TRAINING_OPTIONS_KEY]
    )

    if target_height_m_agl is None:
        scalar_target_names = (
            generator_option_dict[neural_net.SCALAR_TARGET_NAMES_KEY]
        )
        channel_index = scalar_target_names.index(target_field_name)

        actual_values = (
            prediction_dict[prediction_io.SCALAR_TARGETS_KEY][
                example_indices, channel_index
            ]
        )
        predicted_values = (
            prediction_dict[prediction_io.SCALAR_PREDICTIONS_KEY][
                example_indices, channel_index
            ]
        )

        return predicted_values, actual_values

    vector_target_names = (
        generator_option_dict[neural_net.VECTOR_TARGET_NAMES_KEY]
    )
    channel_index = vector_target_names.index(target_field_name)

    height_index = example_io.match_heights(
        heights_m_agl=generator_option_dict[neural_net.HEIGHTS_KEY],
        desired_height_m_agl=target_height_m_agl
    )

    actual_values = (
        prediction_dict[prediction_io.VECTOR_TARGETS_KEY][
            example_indices, height_index, channel_index
        ]
    )
    predicted_values = (
        prediction_dict[prediction_io.VECTOR_PREDICTIONS_KEY][
            example_indices, height_index, channel_index
        ]
    )

    return predicted_values, actual_values


def _plot_saliency_one_example(
        saliency_dict, example_index, model_metadata_dict, legend_suffix,
        output_dir_name):
    """Plots saliency map for one example.

    :param saliency_dict: Dictionary read by `saliency.read_standard_file`.
    :param example_index: Will plot saliency map for example with this array
        index.
    :param model_metadata_dict: Dictionary read by `neural_net.read_metafile`.
    :param legend_suffix: End of figure legend.
    :param output_dir_name: Name of output directory.  Figure will be saved
        here.
    """

    # Housekeeping.
    example_id_string = saliency_dict[saliency.EXAMPLE_IDS_KEY][example_index]
    generator_option_dict = model_metadata_dict[neural_net.TRAINING_OPTIONS_KEY]

    example_dict = {
        example_io.SCALAR_PREDICTOR_NAMES_KEY:
            generator_option_dict[neural_net.SCALAR_PREDICTOR_NAMES_KEY],
        example_io.VECTOR_PREDICTOR_NAMES_KEY:
            generator_option_dict[neural_net.VECTOR_PREDICTOR_NAMES_KEY],
        example_io.HEIGHTS_KEY: generator_option_dict[neural_net.HEIGHTS_KEY],
        example_io.SCALAR_PREDICTOR_VALS_KEY:
            saliency_dict[saliency.SCALAR_SALIENCY_KEY][[example_index], ...],
        example_io.VECTOR_PREDICTOR_VALS_KEY:
            saliency_dict[saliency.VECTOR_SALIENCY_KEY][[example_index], ...]
    }

    scalar_predictor_names = example_dict[example_io.SCALAR_PREDICTOR_NAMES_KEY]
    scalar_saliency_matrix = saliency_dict[saliency.SCALAR_SALIENCY_KEY]

    num_scalar_dim = len(scalar_saliency_matrix.shape) - 1
    num_scalar_predictors = len(scalar_predictor_names)
    legend_string = ''

    if num_scalar_dim == 1:
        for k in range(num_scalar_predictors):
            if k > 0:
                legend_string += '\n'

            legend_string += '{0:s}: {1:.2f}'.format(
                PREDICTOR_NAME_TO_VERBOSE[scalar_predictor_names[k]],
                scalar_saliency_matrix[0, k]
            )

        if legend_suffix != '':
            legend_string += '\n' + legend_suffix
    else:
        legend_string = legend_suffix

    # Plot first set of vector predictors.
    these_flags = numpy.array([
        n in example_dict[example_io.VECTOR_PREDICTOR_NAMES_KEY]
        for n in FIRST_VECTOR_PREDICTOR_NAMES
    ], dtype=bool)

    these_indices = numpy.where(these_flags)[0]

    if len(these_indices) > 0:
        predictor_names = [
            FIRST_VECTOR_PREDICTOR_NAMES[k] for k in these_indices
        ]
        predictor_colours = [
            FIRST_VECTOR_PREDICTOR_COLOURS[k] for k in these_indices
        ]

        handle_dict = profile_plotting.plot_predictors(
            example_dict=example_dict, example_index=0,
            predictor_names=predictor_names,
            predictor_colours=predictor_colours,
            predictor_line_widths=numpy.full(len(these_indices), 2),
            predictor_line_styles=['solid'] * len(these_indices),
            use_log_scale=True, include_units=False, handle_dict=None
        )

        axes_object = handle_dict[profile_plotting.AXES_OBJECTS_KEY][0]

        if legend_string != '':
            axes_object.text(
                0.01, 0.5, legend_string, fontsize=LEGEND_FONT_SIZE, color='k',
                bbox=LEGEND_BOUNDING_BOX_DICT, horizontalalignment='left',
                verticalalignment='center', transform=axes_object.transAxes,
                zorder=1e10
            )

        output_file_name = '{0:s}/{1:s}_first_vector_predictors.jpg'.format(
            output_dir_name, example_id_string.replace('_', '-')
        )
        figure_object = handle_dict[profile_plotting.FIGURE_HANDLE_KEY]

        print('Saving figure to: "{0:s}"...'.format(output_file_name))
        figure_object.savefig(
            output_file_name, dpi=FIGURE_RESOLUTION_DPI, pad_inches=0,
            bbox_inches='tight'
        )
        pyplot.close(figure_object)

    # Plot second set of vector predictors.
    these_flags = numpy.array([
        n in example_dict[example_io.VECTOR_PREDICTOR_NAMES_KEY]
        for n in SECOND_VECTOR_PREDICTOR_NAMES
    ], dtype=bool)

    these_indices = numpy.where(these_flags)[0]

    if len(these_indices) > 0:
        predictor_names = [
            SECOND_VECTOR_PREDICTOR_NAMES[k] for k in these_indices
        ]
        predictor_colours = [
            SECOND_VECTOR_PREDICTOR_COLOURS[k] for k in these_indices
        ]

        handle_dict = profile_plotting.plot_predictors(
            example_dict=example_dict, example_index=0,
            predictor_names=predictor_names,
            predictor_colours=predictor_colours,
            predictor_line_widths=numpy.full(len(these_indices), 2),
            predictor_line_styles=['solid'] * len(these_indices),
            use_log_scale=True, include_units=False, handle_dict=None
        )

        axes_object = handle_dict[profile_plotting.AXES_OBJECTS_KEY][0]

        if legend_string != '':
            axes_object.text(
                0.01, 0.5, legend_string, fontsize=LEGEND_FONT_SIZE, color='k',
                bbox=LEGEND_BOUNDING_BOX_DICT, horizontalalignment='left',
                verticalalignment='center', transform=axes_object.transAxes,
                zorder=1e10
            )

        output_file_name = '{0:s}/{1:s}_second_vector_predictors.jpg'.format(
            output_dir_name, example_id_string.replace('_', '-')
        )
        figure_object = handle_dict[profile_plotting.FIGURE_HANDLE_KEY]

        print('Saving figure to: "{0:s}"...'.format(output_file_name))
        figure_object.savefig(
            output_file_name, dpi=FIGURE_RESOLUTION_DPI, pad_inches=0,
            bbox_inches='tight'
        )
        pyplot.close(figure_object)

    # Plot third set of vector predictors.
    these_flags = numpy.array([
        n in example_dict[example_io.VECTOR_PREDICTOR_NAMES_KEY]
        for n in THIRD_VECTOR_PREDICTOR_NAMES
    ], dtype=bool)

    these_indices = numpy.where(these_flags)[0]

    if len(these_indices) > 0:
        predictor_names = [
            THIRD_VECTOR_PREDICTOR_NAMES[k] for k in these_indices
        ]
        predictor_colours = [
            THIRD_VECTOR_PREDICTOR_COLOURS[k] for k in these_indices
        ]

        handle_dict = profile_plotting.plot_predictors(
            example_dict=example_dict, example_index=0,
            predictor_names=predictor_names,
            predictor_colours=predictor_colours,
            predictor_line_widths=numpy.full(len(these_indices), 2),
            predictor_line_styles=['solid'] * len(these_indices),
            use_log_scale=True, include_units=False, handle_dict=None
        )

        axes_object = handle_dict[profile_plotting.AXES_OBJECTS_KEY][0]

        if legend_string != '':
            axes_object.text(
                0.01, 0.5, legend_string, fontsize=LEGEND_FONT_SIZE, color='k',
                bbox=LEGEND_BOUNDING_BOX_DICT, horizontalalignment='left',
                verticalalignment='center', transform=axes_object.transAxes,
                zorder=1e10
            )

        output_file_name = '{0:s}/{1:s}_third_vector_predictors.jpg'.format(
            output_dir_name, example_id_string.replace('_', '-')
        )
        figure_object = handle_dict[profile_plotting.FIGURE_HANDLE_KEY]

        print('Saving figure to: "{0:s}"...'.format(output_file_name))
        figure_object.savefig(
            output_file_name, dpi=FIGURE_RESOLUTION_DPI, pad_inches=0,
            bbox_inches='tight'
        )
        pyplot.close(figure_object)

    if num_scalar_dim == 1:
        return

    # Plot first set of scalar predictors.
    these_flags = numpy.array([
        n in example_dict[example_io.SCALAR_PREDICTOR_NAMES_KEY]
        for n in FIRST_SCALAR_PREDICTOR_NAMES
    ], dtype=bool)

    these_indices = numpy.where(these_flags)[0]

    if len(these_indices) > 0:
        predictor_names = [
            FIRST_SCALAR_PREDICTOR_NAMES[k] for k in these_indices
        ]
        predictor_colours = [
            FIRST_SCALAR_PREDICTOR_COLOURS[k] for k in these_indices
        ]

        handle_dict = profile_plotting.plot_predictors(
            example_dict=example_dict, example_index=0,
            predictor_names=predictor_names,
            predictor_colours=predictor_colours,
            predictor_line_widths=numpy.full(len(these_indices), 2),
            predictor_line_styles=['solid'] * len(these_indices),
            use_log_scale=True, include_units=False, handle_dict=None
        )

        output_file_name = '{0:s}/{1:s}_first_scalar_predictors.jpg'.format(
            output_dir_name, example_id_string.replace('_', '-')
        )
        figure_object = handle_dict[profile_plotting.FIGURE_HANDLE_KEY]

        print('Saving figure to: "{0:s}"...'.format(output_file_name))
        figure_object.savefig(
            output_file_name, dpi=FIGURE_RESOLUTION_DPI, pad_inches=0,
            bbox_inches='tight'
        )
        pyplot.close(figure_object)

    # Plot second set of scalar predictors.
    these_flags = numpy.array([
        n in example_dict[example_io.SCALAR_PREDICTOR_NAMES_KEY]
        for n in SECOND_SCALAR_PREDICTOR_NAMES
    ], dtype=bool)

    these_indices = numpy.where(these_flags)[0]

    if len(these_indices) > 0:
        predictor_names = [
            SECOND_SCALAR_PREDICTOR_NAMES[k] for k in these_indices
        ]
        predictor_colours = [
            SECOND_SCALAR_PREDICTOR_COLOURS[k] for k in these_indices
        ]

        handle_dict = profile_plotting.plot_predictors(
            example_dict=example_dict, example_index=0,
            predictor_names=predictor_names,
            predictor_colours=predictor_colours,
            predictor_line_widths=numpy.full(len(these_indices), 2),
            predictor_line_styles=['solid'] * len(these_indices),
            use_log_scale=True, include_units=False, handle_dict=None
        )

        output_file_name = '{0:s}/{1:s}_second_scalar_predictors.jpg'.format(
            output_dir_name, example_id_string.replace('_', '-')
        )
        figure_object = handle_dict[profile_plotting.FIGURE_HANDLE_KEY]

        print('Saving figure to: "{0:s}"...'.format(output_file_name))
        figure_object.savefig(
            output_file_name, dpi=FIGURE_RESOLUTION_DPI, pad_inches=0,
            bbox_inches='tight'
        )
        pyplot.close(figure_object)


def _run(saliency_file_name, prediction_file_name, output_dir_name):
    """Plots saliency maps (one for each example).

    This is effectively the main method.

    :param saliency_file_name: See documentation at top of file.
    :param prediction_file_name: Same.
    :param output_dir_name: Same.
    """

    file_system_utils.mkdir_recursive_if_necessary(
        directory_name=output_dir_name
    )

    print('Reading saliency values from: "{0:s}"...'.format(saliency_file_name))
    saliency_dict = saliency.read_standard_file(saliency_file_name)

    target_field_name = saliency_dict[saliency.TARGET_FIELD_KEY]
    target_height_m_agl = saliency_dict[saliency.TARGET_HEIGHT_KEY]
    example_id_strings = saliency_dict[saliency.EXAMPLE_IDS_KEY]
    model_file_name = saliency_dict[saliency.MODEL_FILE_KEY]
    model_metafile_name = neural_net.find_metafile(
        model_dir_name=os.path.split(model_file_name)[0]
    )

    print('Reading model metadata from: "{0:s}"...'.format(model_metafile_name))
    model_metadata_dict = neural_net.read_metafile(model_metafile_name)
    num_examples = len(example_id_strings)

    if target_field_name is None:
        predicted_target_values = [None] * num_examples
        actual_target_values = [None] * num_examples
    else:
        predicted_target_values, actual_target_values = _get_target_values(
            prediction_file_name=prediction_file_name,
            model_metadata_dict=model_metadata_dict,
            example_id_strings=example_id_strings,
            target_field_name=target_field_name,
            target_height_m_agl=target_height_m_agl
        )

    print(SEPARATOR_STRING)

    for i in range(num_examples):
        if target_field_name is None:
            this_legend_suffix = ''
        else:
            this_height_string = (
                '' if target_height_m_agl is None
                else '{0:.2f}-km '.format(METRES_TO_KM * target_height_m_agl)
            )
            this_legend_suffix = (
                'A&P {0:s}{1:s} = {2:.2f}, {3:.2f} {4:s}'
            ).format(
                this_height_string, TARGET_NAME_TO_VERBOSE[target_field_name],
                actual_target_values[i], predicted_target_values[i],
                TARGET_NAME_TO_UNITS[target_field_name]
            )

        _plot_saliency_one_example(
            saliency_dict=saliency_dict, example_index=i,
            model_metadata_dict=model_metadata_dict,
            legend_suffix=this_legend_suffix, output_dir_name=output_dir_name
        )

        if i != num_examples - 1:
            print('\n')


if __name__ == '__main__':
    INPUT_ARG_OBJECT = INPUT_ARG_PARSER.parse_args()

    _run(
        saliency_file_name=getattr(INPUT_ARG_OBJECT, SALIENCY_FILE_ARG_NAME),
        prediction_file_name=getattr(
            INPUT_ARG_OBJECT, PREDICTION_FILE_ARG_NAME
        ),
        output_dir_name=getattr(INPUT_ARG_OBJECT, OUTPUT_DIR_ARG_NAME)
    )
