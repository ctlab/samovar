import argparse
from samovar.parse_annotators import Annotation, match_annotation
import os

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", "-i", type=str, required=True)
parser.add_argument("--output_dir", "-o", type=str, required=False)
parser.add_argument("--true_annotation", "-t", type=str, required=False, default="(?<=taxid:)[0-9]*")
parser.add_argument("--split_sample_name", "-s", type=int, required=False, default=1)
args = parser.parse_args()

def split_sample_name(sample_name: str) -> str:
    return "_".join(sample_name.split("_")[i] for i in range(args.split_sample_name))

# Get sample names
sample_files = set([split_sample_name(gr) for gr in os.listdir(args.input_dir)])

# Process each sample
for sample in sample_files:
    files = [gr for gr in os.listdir(args.input_dir) if split_sample_name(gr) == sample]
    files = [os.path.join(args.input_dir, f) for f in files]

    annotation_dict = {}
    for file in files:
        tool = match_annotation(file)
        if tool is None:
            continue
        annotation_dict[file] = tool

    # Create annotation object
    ann = Annotation(annotation_dict, args.true_annotation)
    ann.correct_level(level="species")

    os.makedirs(args.output_dir, exist_ok=True)
    ann.export(os.path.join(args.output_dir, f"{sample}.annotation.csv"))
    print(f"Exported {sample}")
