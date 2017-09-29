#!/bin/sh

#cp index.php into all subdirectories
find . -mindepth 1 -type d -exec cp index.php {} \;
