import io
from yyc import pipeline
from yyc import scheme
from yyc.utils import data_handle

read_file_path = "../Test_files_for_experimental_validation/Text2.txt"
dna_path = "./output/Text2.dna"
model_path = "./output/Text2.pkl"
write_file_path = "./output/Text2.txt"

if __name__ == "__main__":
    [support_base, rule1, rule2] = ["A", [0, 1, 0, 1], [[1, 1, 0, 0], [1, 0, 0, 1], [1, 1, 0, 0], [1, 1, 0, 0]]]
    tool = scheme.YYC(support_bases=support_base, base_reference=rule1, current_code_matrix=rule2,
                      search_count=100, max_homopolymer=4, max_content=0.6)
    source_stream = io.BytesIO(open(read_file_path, "rb").read())
    dna_stream = io.BytesIO()
    pipeline.encode(
        method=tool,
        input_data=source_stream,
        output_data=dna_stream,
        model_path=model_path,
        need_index=True,
        need_log=True
    )
    with open(dna_path, "wb") as fd:
        dna_stream.seek(0)
        fd.write(dna_stream.read())
    del tool
    dna_stream.seek(0)
    output_stream = io.BytesIO()
    pipeline.decode(
        model_path=model_path,
        input_data=dna_stream,
        output_data=output_stream,
        has_index=True,
        need_log=True
    )
    with open(write_file_path, "wb") as fd:
        output_stream.seek(0)
        fd.write(output_stream.read())

    # compare two file
    matrix_1, _ = data_handle.read_binary_from_all(read_file_path, 120, False)
    matrix_2, _ = data_handle.read_binary_from_all(write_file_path, 120, False)
    print(f"source digital file == target digital file: {matrix_1 == matrix_2}")
