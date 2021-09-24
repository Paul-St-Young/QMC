#!/usr/bin/env python3
import os
#import subprocess as sp
import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt

def main():
  from argparse import ArgumentParser
  parser = ArgumentParser()
  #parser.add_argument('input_filename', type=str)
  parser.add_argument('--verbose', action='store_true')
  args = parser.parse_args()

if __name__ == '__main__':
  main()  # set no global variable
