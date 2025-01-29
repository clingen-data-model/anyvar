#!/usr/bin/env python3
import argparse
import sys

import pysam
import requests


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Load VCF file into AnyVar REST API")
    parser.add_argument("vcf_file", help="VCF file to load")
    parser.add_argument(
        "--anyvar-url", help="Database to load into", default="http://localhost:8000"
    )
    parser.add_argument(
        "--limit", help="Limit number of records to load", type=int, default=None
    )
    return parser.parse_args(argv)


def get_reference_name(variant_file: pysam.VariantFile) -> str:
    reference_headers = [
        record for record in variant_file.header.records if record.key == "reference"
    ]
    return (
        "GRCh38"
        if any("GRCh38" in record.value for record in reference_headers)
        else "GRCh37"
    )


def send_anyvar_put_allele(definition: str, anyvar_endpoint: str):
    body = {"input_type": "Allele", "definition": definition}
    response = requests.put(f"{anyvar_endpoint}/variation", json=body, timeout=10)
    response.raise_for_status()
    return response.json()


def main(argv=sys.argv[1:]):
    args = parse_args(argv)

    with pysam.VariantFile(args.vcf_file, index_filename=None) as f:
        for line_num, variant in enumerate(f):
            if args.limit is not None and line_num >= args.limit:
                break

            print(dir(variant))
            print(f"chrom: {variant.chrom}")
            # print(f"contig: {variant.contig}")
            print(f"pos: {variant.pos}")
            print(f"ref: {variant.ref}")
            print(f"alts: {variant.alts}")
            gnomad_strs = [
                f"{variant.chrom}-{variant.pos}-{variant.ref}-{alt}"
                for alt in variant.alts
            ]
            print(f"gnomad_strs: {gnomad_strs}")
            for gnomad_str in gnomad_strs:
                response = send_anyvar_put_allele(gnomad_str, args.anyvar_url)
                if len(response["messages"]) > 0:
                    raise RuntimeError(f"Error: {response['messages']}")

                vrs_obj = response["object"]
                print(f"vrs_obj: {vrs_obj}")


if __name__ == "__main__":
    main()
