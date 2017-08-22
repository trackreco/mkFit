#!/usr/bin/env python

import os
import argparse
import gpu_setup


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rootdir', dest='rootdir',
                        action='store', default='.',
                        help='Where Makefile.config and Config.h are.')
    args = parser.parse_args()
    if not os.path.isdir(args.rootdir):
        raise RuntimeError('--rootdir should be a directory.')
    return args


def main():
    args = get_args()
    makefile, config_h = gpu_setup.get_file_paths(os.path.abspath(args.rootdir))
    gpu_setup.check_env()
    gpu_setup.backup_original_files([makefile, config_h])

    gpu_setup.setup_files(makefile, config_h)


if __name__ == '__main__':
    main()
