"""
Name: Data Handle

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Function(s):
Conversion of DNA sequences and binary document
"""

import io
import logging
import math
import os
import struct
import typing

from yyc.utils.monitor import Monitor


# noinspection PyProtectedMember
def read_binary_from_all(input_data: typing.Union[str, typing.BinaryIO],
                         segment_length=120, need_log=False):
    """
    introduction: Reading binary matrix from document.

    :param input_data: A stream object or the path of binary file you need to convert.
                       Type: BinaryIO or String.

    :param segment_length: The binary segment length used for DNA sequence generation.
                           Considering current DNA synthesis technique limitation,
                           we usually set 120 as default segment length.

    :param need_log: choose to output log file or not.

    :return matrix: A matrix in which each row represents a binary segment that will be used for DNA sequence generation.
                    Type: two-dimensional list(int)
    """
    if isinstance(input_data, str):
        size = os.path.getsize(input_data)
        logging.debug(f"Read binary matrix from file: {input_data}")
        input_data = open(input_data, mode="rb")
    elif isinstance(input_data, io.BytesIO):
        size = input_data.seek(0, io.SEEK_END)
        input_data.seek(0)

    m = Monitor()

    # Set init storage matrix
    matrix = [[0 for _ in range(segment_length)]
              for _ in range(math.ceil(size * 8 / segment_length))]

    row, col = 0, 0
    for byte_index in range(size):
        if need_log:
            m.output(byte_index + 1, size)
        # Read a file as bytes
        one_byte = input_data.read(1)
        element = list(map(int, "{:08b}".format(ord(one_byte))))
        for item in element:
            matrix[row][col] = item
            col += 1
            if col == segment_length:
                col = 0
                row += 1

    if len(f"{len(matrix):b}") * 7 > segment_length and need_log:
        logging.warn("The proportion of index in whole sequence may be high. \n"
                   "It is recommended to increase the length of output DNA sequences "
                   "or to divide the file into more segment pools")

    input_data.close()
    return matrix, size


# noinspection PyBroadException,PyProtectedMember
def write_all_from_binary(output_data: typing.Union[str, typing.BinaryIO],
                          matrix, size, need_log=False):
    """
    introduction: Writing binary matrix to document.

    :param output_data: A stream object or file path name.
                        Type: BinaryIO or String.

    :param matrix: A matrix in which each row represents a binary segment that will be used for DNA sequence generation.
                   Type: two-dimensional list(int)

    :param size: This refers to file size, to reduce redundant bits when transferring DNA to binary files.
                 Type: int

    :param need_log: choose to output log file or not.
    """
    if isinstance(output_data, str):
        logging.debug(f"Write file from binary matrix: {output_data}")
        output_data = open(output_data, mode="wb+")

    m = Monitor()

    # Change bit to byte (8 -> 1), and write a file as bytes
    bit_index, temp_byte = 0, 0
    for row_id, row in enumerate(matrix):
        if need_log:
            m.output(row_id + 1, len(matrix))
        for cell in row:
            bit_index += 1
            temp_byte *= 2
            temp_byte += cell
            if bit_index == 8:
                if size > 0:
                    output_data.write(struct.pack("B", int(temp_byte)))
                    bit_index = 0
                    temp_byte = 0
                    size -= 1


# noinspection PyBroadException,PyProtectedMember
def read_dna_file(input_data: typing.Union[str, typing.BinaryIO], need_log=False):
    """
    introduction: Reading DNA sequence set from documents.

    :param input_data: A stream object or file path name.
                       Type: BinaryIO or String.

    :return dna_sequences: A corresponding DNA sequence string in which each row acts as a sequence.
                           Type: one-dimensional list(string)

    :param need_log: need output log.
    """
    if isinstance(input_data, str):
        logging.debug(f"Read DNA sequences from file: {input_data}")
        input_data = open(input_data, mode="rb")

    m = Monitor()

    dna_sequences = []

    # Read current file by line
    lines = input_data.readlines()
    for index, line in enumerate(lines):
        if need_log:
            m.output(index + 1, len(lines))
        dna_sequences.append([chr(cell) for cell in line.strip()])

    return dna_sequences


# noinspection PyProtectedMember,PyBroadException
def write_dna_file(output_data: typing.Union[str, typing.BinaryIO], dna_sequences, need_log=False):
    """
    introduction: Writing DNA sequence set to documents.

    :param output_data: A stream object or file path name.
                        Type: BinaryIO or String.

    :param dna_sequences: Generated DNA sequences.
                          Type: one-dimensional list(string)

    :param need_log: choose to output log file or not.
    """
    if isinstance(output_data, str):
        logging.debug(f"Write DNA sequences to binary file: {output_data}")
        output_data = open(output_data, mode="wb")

    m = Monitor()

    for row_id, row in enumerate(dna_sequences):
        if need_log:
            m.output(row_id + 1, len(dna_sequences))
        output_data.write(b"".join(item.encode("utf-8") for item in row) + b"\n")
    return dna_sequences
