#!/usr/bin/env bash

salloc -N 1 --ntasks-per-node=4 --ntasks-per-socket=4 --gres=gpu:1 -t 2:00:00
