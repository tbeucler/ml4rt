"""Methods for building U-nets."""

import keras
from gewittergefahr.gg_utils import error_checking
from gewittergefahr.deep_learning import architecture_utils
from ml4rt.machine_learning import neural_net

NUM_HEIGHTS_KEY = 'num_heights'
NUM_INPUT_CHANNELS_KEY = 'num_input_channels'
INNER_ACTIV_FUNCTION_KEY = 'inner_activ_function_name'
INNER_ACTIV_FUNCTION_ALPHA_KEY = 'inner_activ_function_alpha'
OUTPUT_ACTIV_FUNCTION_KEY = 'output_activ_function_name'
OUTPUT_ACTIV_FUNCTION_ALPHA_KEY = 'output_activ_function_alpha'
L1_WEIGHT_KEY = 'l1_weight'
L2_WEIGHT_KEY = 'l2_weight'
USE_BATCH_NORM_KEY = 'use_batch_normalization'

DEFAULT_ARCHITECTURE_OPTION_DICT = {
    INNER_ACTIV_FUNCTION_KEY: architecture_utils.RELU_FUNCTION_STRING,
    INNER_ACTIV_FUNCTION_ALPHA_KEY: 0.2,
    OUTPUT_ACTIV_FUNCTION_KEY: architecture_utils.RELU_FUNCTION_STRING,
    OUTPUT_ACTIV_FUNCTION_ALPHA_KEY: 0.,
    L1_WEIGHT_KEY: 0.,
    L2_WEIGHT_KEY: 0.001,
    USE_BATCH_NORM_KEY: True
}


def _check_architecture_args(option_dict):
    """Error-checks input arguments for architecture.

    :param option_dict: See doc for `create_model`.
    :return: option_dict: Same as input, except defaults may have been added.
    """

    orig_option_dict = option_dict.copy()
    option_dict = DEFAULT_ARCHITECTURE_OPTION_DICT.copy()
    option_dict.update(orig_option_dict)

    num_heights = option_dict[NUM_HEIGHTS_KEY]
    error_checking.assert_is_integer(num_heights)
    error_checking.assert_is_geq(num_heights, 10)

    num_input_channels = option_dict[NUM_INPUT_CHANNELS_KEY]
    error_checking.assert_is_integer(num_input_channels)
    error_checking.assert_is_geq(num_input_channels, 1)

    l1_weight = option_dict[L1_WEIGHT_KEY]
    error_checking.assert_is_geq(l1_weight, 0.)

    l2_weight = option_dict[L2_WEIGHT_KEY]
    error_checking.assert_is_geq(l2_weight, 0.)

    use_batch_normalization = option_dict[USE_BATCH_NORM_KEY]
    error_checking.assert_is_boolean(use_batch_normalization)

    return option_dict


def make_u_net(option_dict):
    """Creates U-net.

    This method sets up the architecture, loss function, and optimizer -- and
    compiles the model -- but does not train it.

    Architecture taken from:
    https://github.com/zhixuhao/unet/blob/master/model.py

    :param option_dict: Dictionary with the following keys.
    option_dict['num_heights']: Number of height levels.
    option_dict['num_input_channels']: Number of input channels.
    option_dict['inner_activ_function_name']: Name of activation function for
        all inner (non-output) layers.  Must be accepted by
        `architecture_utils.check_activation_function`.
    option_dict['inner_activ_function_alpha']: Alpha (slope parameter) for
        activation function for all inner layers.  Applies only to ReLU and eLU.
    option_dict['output_activ_function_name']: Same as
        `inner_activ_function_name` but for output layers (profiles and
        scalars).
    option_dict['output_activ_function_alpha']: Same as
        `inner_activ_function_alpha` but for output layers (profiles and
        scalars).
    option_dict['l1_weight']: Weight for L_1 regularization.
    option_dict['l2_weight']: Weight for L_2 regularization.
    option_dict['use_batch_normalization']: Boolean flag.  If True, will use
        batch normalization after each inner (non-output) layer.

    :return: model_object: Instance of `keras.models.Model`, with the
        aforementioned architecture.
    """

    # TODO(thunderhoser): Generalize this method a bit.

    option_dict = _check_architecture_args(option_dict=option_dict)

    num_heights = option_dict[NUM_HEIGHTS_KEY]
    num_input_channels = option_dict[NUM_INPUT_CHANNELS_KEY]
    inner_activ_function_name = option_dict[INNER_ACTIV_FUNCTION_KEY]
    inner_activ_function_alpha = option_dict[INNER_ACTIV_FUNCTION_ALPHA_KEY]
    output_activ_function_name = option_dict[OUTPUT_ACTIV_FUNCTION_KEY]
    output_activ_function_alpha = option_dict[OUTPUT_ACTIV_FUNCTION_ALPHA_KEY]
    l1_weight = option_dict[L1_WEIGHT_KEY]
    l2_weight = option_dict[L2_WEIGHT_KEY]
    use_batch_normalization = option_dict[USE_BATCH_NORM_KEY]

    input_layer_object = keras.layers.Input(
        shape=(num_heights, num_input_channels)
    )
    regularizer_object = architecture_utils.get_weight_regularizer(
        l1_weight=l1_weight, l2_weight=l2_weight
    )

    conv_layer1_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = input_layer_object
        else:
            this_input_layer_object = conv_layer1_object

        conv_layer1_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=64,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        conv_layer1_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(conv_layer1_object)

        if use_batch_normalization:
            conv_layer1_object = architecture_utils.get_batch_norm_layer()(
                conv_layer1_object
            )

    pooling_layer1_object = architecture_utils.get_1d_pooling_layer(
        num_rows_in_window=2, num_rows_per_stride=2,
        pooling_type_string=architecture_utils.MAX_POOLING_STRING
    )(conv_layer1_object)

    conv_layer2_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = pooling_layer1_object
        else:
            this_input_layer_object = conv_layer2_object

        conv_layer2_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=128,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        conv_layer2_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(conv_layer2_object)

        if use_batch_normalization:
            conv_layer2_object = architecture_utils.get_batch_norm_layer()(
                conv_layer2_object
            )

    pooling_layer2_object = architecture_utils.get_1d_pooling_layer(
        num_rows_in_window=2, num_rows_per_stride=2,
        pooling_type_string=architecture_utils.MAX_POOLING_STRING
    )(conv_layer2_object)

    conv_layer3_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = pooling_layer2_object
        else:
            this_input_layer_object = conv_layer3_object

        conv_layer3_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=256,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        conv_layer3_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(conv_layer3_object)

        if use_batch_normalization:
            conv_layer3_object = architecture_utils.get_batch_norm_layer()(
                conv_layer3_object
            )

    pooling_layer3_object = architecture_utils.get_1d_pooling_layer(
        num_rows_in_window=2, num_rows_per_stride=2,
        pooling_type_string=architecture_utils.MAX_POOLING_STRING
    )(conv_layer3_object)

    conv_layer4_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = pooling_layer3_object
        else:
            this_input_layer_object = conv_layer4_object

        conv_layer4_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=512,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        conv_layer4_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(conv_layer4_object)

        conv_layer4_object = architecture_utils.get_dropout_layer(
            dropout_fraction=0.5
        )(conv_layer4_object)

        if use_batch_normalization:
            conv_layer4_object = architecture_utils.get_batch_norm_layer()(
                conv_layer4_object
            )

    pooling_layer4_object = architecture_utils.get_1d_pooling_layer(
        num_rows_in_window=2, num_rows_per_stride=2,
        pooling_type_string=architecture_utils.MAX_POOLING_STRING
    )(conv_layer4_object)

    conv_layer5_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = pooling_layer4_object
        else:
            this_input_layer_object = conv_layer5_object

        conv_layer5_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=1024,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        conv_layer5_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(conv_layer5_object)

        conv_layer5_object = architecture_utils.get_dropout_layer(
            dropout_fraction=0.5
        )(conv_layer5_object)

        if use_batch_normalization:
            conv_layer5_object = architecture_utils.get_batch_norm_layer()(
                conv_layer5_object
            )

    this_layer_object = keras.layers.UpSampling1D(size=2)

    upconv_layer4_object = architecture_utils.get_1d_conv_layer(
        num_kernel_rows=2, num_rows_per_stride=1, num_filters=512,
        padding_type_string=architecture_utils.YES_PADDING_STRING,
        weight_regularizer=regularizer_object
    )(this_layer_object(conv_layer5_object))

    merged_layer4_object = keras.layers.Concatenate(axis=-1)(
        [conv_layer4_object, upconv_layer4_object]
    )

    second_conv_layer4_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = merged_layer4_object
        else:
            this_input_layer_object = second_conv_layer4_object

        second_conv_layer4_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=512,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        second_conv_layer4_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(second_conv_layer4_object)

        if use_batch_normalization:
            second_conv_layer4_object = (
                architecture_utils.get_batch_norm_layer()(
                    second_conv_layer4_object
                )
            )

    this_layer_object = keras.layers.UpSampling1D(size=2)

    upconv_layer3_object = architecture_utils.get_1d_conv_layer(
        num_kernel_rows=2, num_rows_per_stride=1, num_filters=256,
        padding_type_string=architecture_utils.YES_PADDING_STRING,
        weight_regularizer=regularizer_object
    )(this_layer_object(second_conv_layer4_object))

    merged_layer3_object = keras.layers.Concatenate(axis=-1)(
        [conv_layer3_object, upconv_layer3_object]
    )

    second_conv_layer3_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = merged_layer3_object
        else:
            this_input_layer_object = second_conv_layer3_object

        second_conv_layer3_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=256,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        second_conv_layer3_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(second_conv_layer3_object)

        if use_batch_normalization:
            second_conv_layer3_object = (
                architecture_utils.get_batch_norm_layer()(
                    second_conv_layer3_object
                )
            )

    this_layer_object = keras.layers.UpSampling1D(size=2)

    upconv_layer2_object = architecture_utils.get_1d_conv_layer(
        num_kernel_rows=2, num_rows_per_stride=1, num_filters=128,
        padding_type_string=architecture_utils.YES_PADDING_STRING,
        weight_regularizer=regularizer_object
    )(this_layer_object(second_conv_layer3_object))

    merged_layer2_object = keras.layers.Concatenate(axis=-1)(
        [conv_layer2_object, upconv_layer2_object]
    )

    second_conv_layer2_object = None

    for i in range(2):
        if i == 0:
            this_input_layer_object = merged_layer2_object
        else:
            this_input_layer_object = second_conv_layer2_object

        second_conv_layer2_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1, num_filters=128,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        second_conv_layer2_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(second_conv_layer2_object)

        if use_batch_normalization:
            second_conv_layer2_object = (
                architecture_utils.get_batch_norm_layer()(
                    second_conv_layer2_object
                )
            )

    this_layer_object = keras.layers.UpSampling1D(size=2)

    upconv_layer1_object = architecture_utils.get_1d_conv_layer(
        num_kernel_rows=2, num_rows_per_stride=1, num_filters=64,
        padding_type_string=architecture_utils.YES_PADDING_STRING,
        weight_regularizer=regularizer_object
    )(this_layer_object(second_conv_layer2_object))

    merged_layer1_object = keras.layers.Concatenate(axis=-1)(
        [conv_layer1_object, upconv_layer1_object]
    )

    second_conv_layer1_object = None

    for i in range(3):
        if i == 0:
            this_input_layer_object = merged_layer1_object
        else:
            this_input_layer_object = second_conv_layer1_object

        second_conv_layer1_object = architecture_utils.get_1d_conv_layer(
            num_kernel_rows=3, num_rows_per_stride=1,
            num_filters=6 if i == 2 else 3,
            padding_type_string=architecture_utils.YES_PADDING_STRING,
            weight_regularizer=regularizer_object
        )(this_input_layer_object)

        second_conv_layer1_object = architecture_utils.get_activation_layer(
            activation_function_string=inner_activ_function_name,
            alpha_for_relu=inner_activ_function_alpha,
            alpha_for_elu=inner_activ_function_alpha
        )(second_conv_layer1_object)

        if use_batch_normalization:
            second_conv_layer1_object = (
                architecture_utils.get_batch_norm_layer()(
                    second_conv_layer1_object
                )
            )

    second_conv_layer1_object = architecture_utils.get_1d_conv_layer(
        num_kernel_rows=1, num_rows_per_stride=1, num_filters=3,
        padding_type_string=architecture_utils.YES_PADDING_STRING,
        weight_regularizer=regularizer_object
    )(second_conv_layer1_object)

    second_conv_layer1_object = architecture_utils.get_activation_layer(
        activation_function_string=output_activ_function_name,
        alpha_for_relu=output_activ_function_alpha,
        alpha_for_elu=output_activ_function_alpha
    )(second_conv_layer1_object)

    model_object = keras.models.Model(
        input=input_layer_object, output=second_conv_layer1_object
    )

    model_object.compile(
        loss=keras.losses.mse, optimizer=keras.optimizers.Adam(),
        metrics=neural_net.METRIC_FUNCTION_LIST
    )

    model_object.summary()
    return model_object