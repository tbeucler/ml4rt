"""Averages many examples."""

import argparse
import numpy
from ml4rt.io import example_io
from ml4rt.utils import misc as misc_utils
from ml4rt.scripts import make_saliency_maps

TIME_FORMAT = '%Y-%m-%d-%H%M%S'

EXAMPLE_FILE_ARG_NAME = make_saliency_maps.EXAMPLE_FILE_ARG_NAME
NUM_EXAMPLES_ARG_NAME = make_saliency_maps.NUM_EXAMPLES_ARG_NAME
EXAMPLE_DIR_ARG_NAME = make_saliency_maps.EXAMPLE_DIR_ARG_NAME
EXAMPLE_ID_FILE_ARG_NAME = make_saliency_maps.EXAMPLE_ID_FILE_ARG_NAME
USE_PMM_ARG_NAME = 'use_pmm'
MAX_PERCENTILE_ARG_NAME = 'max_pmm_percentile_level'
OUTPUT_FILE_ARG_NAME = 'output_file_name'

EXAMPLE_FILE_HELP_STRING = make_saliency_maps.EXAMPLE_FILE_HELP_STRING
NUM_EXAMPLES_HELP_STRING = make_saliency_maps.NUM_EXAMPLES_HELP_STRING
EXAMPLE_DIR_HELP_STRING = make_saliency_maps.EXAMPLE_DIR_HELP_STRING
EXAMPLE_ID_FILE_HELP_STRING = make_saliency_maps.EXAMPLE_ID_FILE_HELP_STRING
USE_PMM_HELP_STRING = (
    'Boolean flag.  If 1 (0), will use probability-matched (arithmetic) means '
    'for vertical profiles.'
)
MAX_PERCENTILE_HELP_STRING = (
    '[used only if `{0:s}` = 1] Max percentile level for probability-matched '
    'means.'
)
OUTPUT_FILE_HELP_STRING = (
    'Path to output file.  Average of all examples will be written here by '
    '`example_io.write_file`.'
)

INPUT_ARG_PARSER = argparse.ArgumentParser()
INPUT_ARG_PARSER.add_argument(
    '--' + EXAMPLE_FILE_ARG_NAME, type=str, required=False, default='',
    help=EXAMPLE_FILE_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + NUM_EXAMPLES_ARG_NAME, type=int, required=False, default=-1,
    help=NUM_EXAMPLES_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + EXAMPLE_DIR_ARG_NAME, type=str, required=False, default='',
    help=EXAMPLE_DIR_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + EXAMPLE_ID_FILE_ARG_NAME, type=str, required=False, default='',
    help=EXAMPLE_ID_FILE_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + USE_PMM_ARG_NAME, type=int, required=True, help=USE_PMM_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + MAX_PERCENTILE_ARG_NAME, type=float, required=False, default=99.,
    help=MAX_PERCENTILE_HELP_STRING
)
INPUT_ARG_PARSER.add_argument(
    '--' + OUTPUT_FILE_ARG_NAME, type=str, required=True,
    help=OUTPUT_FILE_HELP_STRING
)


def _run(example_file_name, num_examples, example_dir_name,
         example_id_file_name, use_pmm, max_pmm_percentile_level,
         output_file_name):
    """Averages many examples.

    This is effectively the main method.

    :param example_file_name: Same.
    :param num_examples: Same.
    :param example_dir_name: Same.
    :param example_id_file_name: Same.
    :param use_pmm: Same.
    :param max_pmm_percentile_level: Same.
    :param output_file_name: Same.
    """

    use_specific_ids = example_file_name == ''

    if use_specific_ids:
        print('Reading desired example IDs from: "{0:s}"...'.format(
            example_id_file_name
        ))
        example_id_strings = (
            misc_utils.read_example_ids_from_netcdf(example_id_file_name)
        )

        valid_times_unix_sec = example_io.parse_example_ids(example_id_strings)[
            example_io.VALID_TIMES_KEY
        ]
        example_file_names = example_io.find_many_files(
            example_dir_name=example_dir_name,
            first_time_unix_sec=numpy.min(valid_times_unix_sec),
            last_time_unix_sec=numpy.max(valid_times_unix_sec)
        )

        num_files = len(example_file_names)
        example_dicts = [dict()] * num_files

        for i in range(num_files):
            print('Reading data from: "{0:s}"...'.format(example_file_names[i]))
            example_dicts[i] = example_io.read_file(example_file_names[i])

        example_dict = example_io.concat_examples(example_dicts)

        all_id_strings = example_io.create_example_ids(example_dict)
        good_indices = example_io.find_examples(
            all_id_strings=all_id_strings,
            desired_id_strings=example_id_strings, allow_missing=False
        )
        example_dict = example_io.subset_by_index(
            example_dict=example_dict, desired_indices=good_indices
        )
    else:
        print('Reading data from: "{0:s}"...'.format(example_file_name))
        example_dict = example_io.read_file(example_file_name)
        example_dict = example_io.reduce_sample_size(
            example_dict=example_dict, num_examples_to_keep=num_examples
        )
    
    print(example_dict[example_io.VALID_TIMES_KEY])
    num_examples = len(
        example_dict[example_io.VALID_TIMES_KEY]
    )

    print('Averaging {0:d} examples...'.format(num_examples))
    mean_example_dict = example_io.average_examples(
        example_dict=example_dict, use_pmm=use_pmm,
        max_pmm_percentile_level=max_pmm_percentile_level
    )

    print('Writing mean example to: "{0:s}"...'.format(output_file_name))
    example_io.write_file(
        example_dict=mean_example_dict, netcdf_file_name=output_file_name
    )


if __name__ == '__main__':
    INPUT_ARG_OBJECT = INPUT_ARG_PARSER.parse_args()

    _run(
        example_file_name=getattr(INPUT_ARG_OBJECT, EXAMPLE_FILE_ARG_NAME),
        num_examples=getattr(INPUT_ARG_OBJECT, NUM_EXAMPLES_ARG_NAME),
        example_dir_name=getattr(INPUT_ARG_OBJECT, EXAMPLE_DIR_ARG_NAME),
        example_id_file_name=getattr(
            INPUT_ARG_OBJECT, EXAMPLE_ID_FILE_ARG_NAME
        ),
        use_pmm=bool(getattr(INPUT_ARG_OBJECT, USE_PMM_ARG_NAME)),
        max_pmm_percentile_level=getattr(
            INPUT_ARG_OBJECT, MAX_PERCENTILE_ARG_NAME
        ),
        output_file_name=getattr(INPUT_ARG_OBJECT, OUTPUT_FILE_ARG_NAME)
    )
