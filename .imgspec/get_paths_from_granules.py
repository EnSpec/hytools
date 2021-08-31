#! /usr/bin/python

import argparse
import glob
import os
import sys
import tarfile


def parse_args():
    parser = argparse.ArgumentParser()
    products = ["obs_ort", "rfl", "topo_coeffs", "brdf_coeffs", "topo_brdf"]
    formats = ["envi"]
    parser.add_argument("-p", "--product",
                        help=("Choose one of the following product types: " + ", ".join(products)))
    parser.add_argument("-f", "--format",
                        help=("Choose one of the following formats: " + ", ".join(formats)))
    args = parser.parse_args()

    if args.product:
        if args.product not in products:
            print("ERROR: Product \"%s\" is not a valid product choice." % args.product)
            sys.exit(1)
    if args.format:
        if args.format not in formats:
            print("ERROR: Format \"%s\" is not a valid format choice." % f)
            sys.exit(1)
    return args


def main():
    args = parse_args()

    # Unzip and untar granules
    input_dir = "input"
    granule_paths = glob.glob(os.path.join(input_dir, "*.tar.gz"))
    for g in granule_paths:
        tar_file = tarfile.open(g)
        tar_file.extractall(input_dir)
        tar_file.close()
        os.remove(g)

    # Get paths based on product type file matching
    # TODO: Add support for multiple formats
    if args.product == "obs_ort":
        paths = glob.glob(os.path.join(input_dir, "*rdn*", "*obs_ort"))
    elif args.product == "rfl":
        paths = glob.glob(os.path.join(input_dir, "*", "*rfl*img"))
        paths += glob.glob(os.path.join(input_dir, "*", "*corr*img"))
    elif args.product == "topo_coeffs":
        paths = glob.glob(os.path.join(input_dir, "*topo_coeffs*", "*topo_coeffs*json"))
    elif args.product == "brdf_coeffs":
        paths = glob.glob(os.path.join(input_dir, "*brdf_coeffs*", "*brdf_coeffs*json"))
    elif args.product == "topo_brdf":
        paths = glob.glob(os.path.join(input_dir, "*", "*topo_brdf"))
    print(",".join(paths))


if __name__ == "__main__":
    main()
