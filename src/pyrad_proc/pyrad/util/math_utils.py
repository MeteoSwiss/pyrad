import math


def bit_pack(values, max_value, num_values):
    """
    Pack a list of integers into a single integer using bit manipulation.

    Parameters
    ----------
    values : list of int
        The values to pack. Each value must be <= `max_value`.
        If fewer than `num_values` are given, values are right-padded with 0s.
        If more than `num_values` are given, a ValueError is raised.
    max_value : int
        The maximum possible value any item in `values` can take.
        Determines how many bits are allocated per value.
    num_values : int
        The exact number of values to pack. Controls total bit length.

    Returns
    -------
    int
        A single integer encoding all values.

    Raises
    ------
    ValueError
        If any value in `values` exceeds `max_value`, or if more than `num_values`
        are provided.

    """
    bits = math.ceil(math.log2(max_value + 1))
    if len(values) > num_values:
        raise ValueError("Too many values provided")
    padded_values = values + [0] * (num_values - len(values))

    packed = 0
    for v in padded_values:
        if v > max_value:
            raise ValueError(f"Value {v} exceeds max allowed {max_value}")
        packed = (packed << bits) | v
    return packed


def bit_unpack(packed, max_value, num_values):
    """
    Unpack a single integer into a list of integers using bit manipulation.

    Parameters
    ----------
    packed : int
        The packed integer to decode.
    max_value : int
        The maximum possible value any original item could have taken.
        Determines how many bits were allocated per value.
    num_values : int
        The number of values to extract.

    Returns
    -------
    list of int
        The unpacked values in their original order.

    """
    bits = math.ceil(math.log2(max_value + 1))
    mask = (1 << bits) - 1
    values = []
    for _ in range(num_values):
        values.append(packed & mask)
        packed >>= bits
    return list(reversed(values))
