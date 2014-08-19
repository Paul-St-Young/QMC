#!/usr/bin/env python

import sys

from gamess_analyzer import GamessAnalyzer
ana=GamessAnalyzer(sys.argv[1])
ana.analyze()
